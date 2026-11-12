"""Concatenate GEddySoft daily HDF5 outputs into a single NetCDF file.

This script reads daily GEddySoft output files stored as HDF5 and builds one
time-continuous xarray Dataset, then writes it to NetCDF3 (scipy backend) for
maximum compatibility.

Configuration ("User inputs") is done directly in the block at the end of this
file under:
    if __name__ == '__main__':

Expected HDF5 structure
-----------------------
- The file contains a top-level group either:
  - TRACER/<i>/... (PTR-TOF-MS style), or
  - IRGA/<i>/... (QCL/GHG style)
  where <i> are numeric channel keys ('0', '1', ...).
- Variables that are per-channel (e.g. flux) are exported with dimensions:
  - (time, mz) for TRACER, or
  - (time, irga) for IRGA.

Time window selection
---------------------
- start_time/end_time are day strings 'YYYY-MM-DD'.
- end_time is inclusive (the script selects until the next day, exclusive).
- Files outside the requested day window are not opened (based on yyyymmdd
  found in filenames), and the final Dataset is also filtered by time.

Optional external MET CSV
-------------------------
- If enabled, an external MET CSV is read and merged on the 'time' coordinate.
- Timestamps are floored to seconds and duplicate timestamps are removed to
  avoid xarray alignment errors.

Output
------
- The output NetCDF is written as NetCDF3 64-bit:
    ds.to_netcdf(..., engine='scipy', format='NETCDF3_64BIT')
"""

import re
from pathlib import Path

import h5py
import numpy as np
import xarray as xr
import pandas as pd


def _split_ipt_like(values: np.ndarray, description: str) -> tuple[list[str], list[np.ndarray]]:
    match = re.search(r"\((.*?)\)", str(description))
    if not match:
        return [], []
    names = [item.strip() for item in match.group(1).split(',')]
    if not names:
        return [], []
    return names, [values[:, i] for i in range(min(values.shape[1], len(names)))]


def _irga_name_to_coord(name_value) -> str:
    if isinstance(name_value, (bytes, bytearray)):
        s = name_value.decode(errors='replace')
    else:
        s = str(name_value)

    if s == 'H2O dry':
        return '18'
    if s == 'N2O dry':
        return '44'
    return s


def _decode_time_to_datetime64(time_values: np.ndarray) -> np.ndarray:
    if np.issubdtype(time_values.dtype, np.bytes_):
        time_values = [t.decode() for t in time_values]
    elif np.issubdtype(time_values.dtype, np.str_):
        time_values = list(time_values)
    else:
        time_values = time_values.astype(str).tolist()

    time_cleaned = [t[:10] + ' ' + t[11:].replace('-', ':') for t in time_values]
    return np.array(time_cleaned, dtype='datetime64[s]')


def _decode_mz_to_float(mz_values: np.ndarray) -> np.ndarray:
    if mz_values.dtype.kind in {'S', 'O', 'U'}:
        return np.array([
            float(m.decode()) if isinstance(m, (bytes, bytearray)) else float(m) for m in mz_values
        ], dtype=float)
    return mz_values.astype(float)


def _decode_mz_to_float_irga_safe(mz_values: np.ndarray, *, filepath: Path) -> np.ndarray:
    try:
        return _decode_mz_to_float(mz_values)
    except (TypeError, ValueError):
        bad: list[str] = []
        for v in mz_values.tolist():
            try:
                float(v.decode() if isinstance(v, (bytes, bytearray)) else v)
            except (TypeError, ValueError):
                bad.append(str(v.decode(errors='replace') if isinstance(v, (bytes, bytearray)) else v))

        bad_unique = sorted(set(bad))
        raise ValueError(
            "IRGA compound names could not be converted to numeric mz values. "
            "This breaks .sel(mz=..., method='nearest'). "
            f"File: {filepath}. Non-numeric values: {bad_unique}"
        )


def _list_daily_files(folder: Path, prefix: str, start: str | None, end: str | None) -> list[Path]:
    return sorted(folder.glob(f"{prefix}*.hdf5"))


def _filter_files_by_time_window(files: list[Path], start_time: str | None, end_time: str | None) -> list[Path]:
    if start_time is None and end_time is None:
        return files

    start_day = np.datetime64(start_time) if start_time is not None and start_time != '' else None
    end_day = np.datetime64(end_time) if end_time is not None and end_time != '' else None

    date_re = re.compile(r"(\d{8})")
    kept: list[Path] = []

    for p in files:
        m = date_re.search(p.name)
        if not m:
            kept.append(p)
            continue

        day_str = m.group(1)
        day = np.datetime64(f"{day_str[0:4]}-{day_str[4:6]}-{day_str[6:8]}")

        if start_day is not None and day < start_day:
            continue
        if end_day is not None and day > end_day:
            continue
        kept.append(p)

    return kept


def _read_single_file(filepath: Path, *, mode: str) -> xr.Dataset:
    with h5py.File(filepath, mode='r') as results:
        time = _decode_time_to_datetime64(results['time'][()])
        freq = results['freq'][()]

        if mode == 'auto':
            if 'TRACER' in results:
                mode = 'tracer'
            elif 'IRGA' in results:
                mode = 'irga'
            else:
                raise ValueError(f"Missing TRACER/IRGA groups in {filepath}")

        if mode == 'tracer':
            channel_group_name = 'TRACER'
            channel_keys = [k for k in results[channel_group_name].keys() if str(k).isdigit()]
            if not channel_keys:
                raise ValueError(f"No numeric TRACER groups found in {filepath}")
            channel_coord = _decode_mz_to_float(np.array([results[channel_group_name][k]['mz'][()] for k in channel_keys]))
        elif mode == 'irga':
            channel_group_name = 'IRGA'
            channel_keys = [k for k in results[channel_group_name].keys() if str(k).isdigit()]
            if not channel_keys:
                raise ValueError(f"No numeric IRGA groups found in {filepath}")
            channel_coord_raw = np.array([
                _irga_name_to_coord(results[channel_group_name][k]['name'][()]) for k in channel_keys
            ], dtype=object)
            channel_coord = _decode_mz_to_float_irga_safe(channel_coord_raw, filepath=filepath)
        else:
            raise ValueError(f"Unsupported mode: {mode}")

        channel_time_vars: dict[str, list[np.ndarray]] = {}
        channel_freq_vars: dict[str, list[np.ndarray]] = {}
        channel_qaqc_time_vars: dict[str, list[np.ndarray]] = {}
        channel_root_time_vars: dict[str, np.ndarray] = {}

        for key in channel_keys:
            g = results[channel_group_name][key]

            for var in g.keys():
                if var in {'mz', 'name', 'qaqc'}:
                    continue
                values = g[var][()]
                if values.ndim == 2 and values.shape == (len(time), len(freq)):
                    channel_freq_vars.setdefault(var, []).append(values)
                elif values.ndim == 1 and values.shape == (len(time),):
                    channel_time_vars.setdefault(var, []).append(values)

            if 'qaqc' in g:
                qaqc = g['qaqc']
                for var in qaqc.keys():
                    if var.startswith('IPT'):
                        values = qaqc[var][()]
                        if values.ndim == 2 and values.shape[0] == len(time):
                            description = qaqc[var].attrs.get('description', '')
                            names, cols = _split_ipt_like(values, description)
                            if names:
                                for name, col in zip(names, cols):
                                    channel_qaqc_time_vars.setdefault(f"{var}_{name}", []).append(col)
                            else:
                                channel_qaqc_time_vars.setdefault(var, []).append(values)
                        continue

                    values = qaqc[var][()]
                    if values.ndim == 1 and values.shape == (len(time),):
                        channel_qaqc_time_vars.setdefault(var, []).append(values)

        for var in results[channel_group_name].keys():
            if str(var).isdigit() or var in {'qaqc'}:
                continue
            obj = results[channel_group_name][var]
            if isinstance(obj, h5py.Dataset):
                values = obj[()]
                if values.ndim == 1 and values.shape == (len(time),):
                    channel_root_time_vars[var] = values

        met_time_vars: dict[str, np.ndarray] = {}
        met_freq_vars: dict[str, np.ndarray] = {}
        met_qaqc_time_vars: dict[str, np.ndarray] = {}

        if 'MET' in results:
            met_group = results['MET']
            for var in met_group.keys():
                if var == 'qaqc':
                    continue
                obj = met_group[var]
                if not isinstance(obj, h5py.Dataset):
                    continue
                values = obj[()]
                if values.ndim == 2 and values.shape == (len(time), len(freq)):
                    met_freq_vars[var] = values
                elif values.ndim == 1 and values.shape == (len(time),):
                    met_time_vars[var] = values

            if 'qaqc' in met_group:
                qaqc = met_group['qaqc']
                for var in qaqc.keys():
                    if var.startswith('IPT'):
                        values = qaqc[var][()]
                        if values.ndim == 2 and values.shape[0] == len(time):
                            description = qaqc[var].attrs.get('description', '')
                            names, cols = _split_ipt_like(values, description)
                            if names:
                                for name, col in zip(names, cols):
                                    met_qaqc_time_vars[f"{var}_{name}"] = col
                            else:
                                met_qaqc_time_vars[var] = values
                        continue

                    obj = qaqc[var]
                    if not isinstance(obj, h5py.Dataset):
                        continue
                    values = obj[()]
                    if values.ndim == 1 and values.shape == (len(time),):
                        met_qaqc_time_vars[var] = values

        coord_name = 'mz' if mode == 'tracer' else 'irga'

        data_vars_xr: dict[str, tuple] = {}

        for var, per_chan in channel_time_vars.items():
            data_np = np.stack(per_chan, axis=0)  # (chan,time)
            attrs = {
                'name': results[channel_group_name][channel_keys[0]][var].attrs.get('description', ''),
                'units': results[channel_group_name][channel_keys[0]][var].attrs.get('units', ''),
            }
            data_vars_xr[var] = (("time", coord_name), data_np.T, attrs)

        for var, per_chan in channel_freq_vars.items():
            data_np = np.stack(per_chan, axis=0)  # (chan,time,freq)
            attrs = {
                'name': results[channel_group_name][channel_keys[0]][var].attrs.get('description', ''),
                'units': results[channel_group_name][channel_keys[0]][var].attrs.get('units', ''),
            }
            data_vars_xr[var] = (("time", coord_name, "freq"), data_np.transpose(1, 0, 2), attrs)

        for var, per_chan in channel_qaqc_time_vars.items():
            data_np = np.stack(per_chan, axis=0)  # (chan,time)
            attrs = {'name': var, 'units': ''}
            if var.startswith('IPT'):
                attrs['units'] = '%'
            else:
                base = var.split('_')[0]
                if 'qaqc' in results[channel_group_name][channel_keys[0]] and base in results[channel_group_name][channel_keys[0]]['qaqc']:
                    attrs = {
                        'name': results[channel_group_name][channel_keys[0]]['qaqc'][base].attrs.get('description', var),
                        'units': results[channel_group_name][channel_keys[0]]['qaqc'][base].attrs.get('units', ''),
                    }
            data_vars_xr[var] = (("time", coord_name), data_np.T, attrs)

        for var, values in channel_root_time_vars.items():
            attrs = {
                'name': results[channel_group_name][var].attrs.get('description', ''),
                'units': results[channel_group_name][var].attrs.get('units', ''),
            }
            data_vars_xr[var] = (("time",), values, attrs)

        for var, values in met_time_vars.items():
            attrs = {
                'name': results['MET'][var].attrs.get('description', ''),
                'units': results['MET'][var].attrs.get('units', ''),
            }
            data_vars_xr[var] = (("time",), values, attrs)

        for var, values in met_freq_vars.items():
            attrs = {
                'name': results['MET'][var].attrs.get('description', ''),
                'units': results['MET'][var].attrs.get('units', ''),
            }
            data_vars_xr[var] = (("time", "freq"), values, attrs)

        for var, values in met_qaqc_time_vars.items():
            attrs = {'name': var, 'units': ''}
            if var.startswith('IPT'):
                attrs['units'] = '%'
            else:
                base = var.split('_')[0]
                if 'MET' in results and 'qaqc' in results['MET'] and base in results['MET']['qaqc']:
                    attrs = {
                        'name': results['MET']['qaqc'][base].attrs.get('description', var),
                        'units': results['MET']['qaqc'][base].attrs.get('units', ''),
                    }
            data_vars_xr[var] = (("time",), values, attrs)

        coords = {"time": time, "freq": freq, coord_name: channel_coord}
        ds = xr.Dataset(data_vars_xr, coords=coords)
        if mode == 'tracer':
            ds = ds.sortby('mz')
        elif mode == 'irga':
            ds = ds.rename({'irga': 'mz'})
        return ds


def _read_external_met_csv(met_csv_path: str) -> xr.Dataset:
    df = pd.read_csv(met_csv_path)
    if df.empty:
        raise ValueError(f"External MET file is empty: {met_csv_path}")

    if 'TIMESTAMP_END' in df.columns:
        time_col = 'TIMESTAMP_END'
    elif 'TIMESTAMP' in df.columns:
        time_col = 'TIMESTAMP'
    else:
        time_col = df.columns[0]

    times = pd.to_datetime(df[time_col], errors='coerce')
    df = df.loc[~times.isna()].copy()
    times = times.loc[~times.isna()]

    # Align external MET timestamps to the same resolution as GEddySoft outputs (seconds).
    # This must happen BEFORE dropping duplicates, otherwise distinct ns timestamps can
    # collapse to identical seconds and break xarray alignment.
    times = times.dt.floor('s')

    df = df.drop(columns=[time_col])
    df.index = times
    df = df[~df.index.duplicated(keep='first')]

    time64 = df.index.to_numpy(dtype='datetime64[s]')

    data_vars = {}
    for col in df.columns:
        if pd.api.types.is_numeric_dtype(df[col]):
            values = df[col].to_numpy()
        else:
            values = df[col].astype(str).to_numpy()
        data_vars[col] = (("time",), values)

    return xr.Dataset(data_vars, coords={"time": time64})


def _drop_duplicate_times(ds: xr.Dataset) -> xr.Dataset:
    if 'time' not in ds.coords:
        return ds
    time_values = ds['time'].values
    if time_values.size == 0:
        return ds

    # Ensure monotonic increasing time, then drop duplicates (keep first occurrence)
    order = np.argsort(time_values)
    ds_sorted = ds.isel(time=order)
    time_sorted = ds_sorted['time'].values

    _, unique_indices = np.unique(time_sorted, return_index=True)
    if unique_indices.size == time_sorted.size:
        return ds_sorted
    return ds_sorted.isel(time=np.sort(unique_indices))


def _count_duplicate_times(ds: xr.Dataset) -> int:
    if 'time' not in ds.coords:
        return 0
    t = ds['time'].values
    if t.size == 0:
        return 0
    _, counts = np.unique(t, return_counts=True)
    return int(np.sum(counts > 1))


def _apply_time_window(ds: xr.Dataset, start_time: str | None, end_time: str | None) -> xr.Dataset:
    if start_time is None and end_time is None:
        return ds

    # Day-level selection:
    # - start_time/end_time are expected as 'YYYY-MM-DD'
    # - end day is inclusive (we slice until next day, exclusive)
    start = np.datetime64(start_time) if start_time is not None and start_time != '' else None
    end = np.datetime64(end_time) if end_time is not None and end_time != '' else None
    end_exclusive = (end + np.timedelta64(1, 'D')) if end is not None else None

    if start is not None and end_exclusive is not None:
        return ds.sel(time=slice(start, end_exclusive))
    if start is not None:
        return ds.sel(time=slice(start, None))
    return ds.sel(time=slice(None, end_exclusive))


def concat_geddysoft_hdf5(
    results_folder: str,
    prefix: str,
    output_nc: str,
    *,
    mode: str = 'auto',
    start_time: str | None = None,
    end_time: str | None = None,
) -> xr.Dataset:
    folder = Path(results_folder)
    if not folder.is_dir():
        raise FileNotFoundError(f"results_folder not found: {results_folder}")

    files = _list_daily_files(folder, prefix, None, None)
    files = _filter_files_by_time_window(files, start_time, end_time)
    if not files:
        raise FileNotFoundError(f"No HDF5 files found for prefix={prefix} in {results_folder}")

    datasets: list[xr.Dataset] = []
    for p in files:
        print(f"Read file: {p.name}")
        datasets.append(_read_single_file(p, mode=mode))

    if not datasets:
        raise RuntimeError("No files were readable after size filtering")

    ds = xr.concat(datasets, dim='time')
    ds = _drop_duplicate_times(ds)
    if 'mz' in ds.coords:
        ds = ds.sortby('mz')

    return ds


if __name__ == '__main__':
    base_dir = Path(__file__).resolve().parent

    # User inputs
    # -----------
    # HDF5 daily outputs folder (relative to this script location)
    results_folder = str((base_dir / '..' / 'output').resolve())

    # Daily results filename prefix. Expected files look like:
    #   <prefix><yyyymmdd>.hdf5
    prefix = 'QCL_BE-Lon_fluxes_'

    # NetCDF output path
    output_nc = str((base_dir / '..' / 'output' / 'QCL_BE-Lon_fluxes_concat.nc').resolve())

    # Which group to read from the daily HDF5 files:
    # - 'auto'  : try TRACER first, else IRGA
    # - 'tracer': read TRACER group (PTR-TOF-MS style)
    # - 'irga'  : read IRGA group (GHG/QCL style)
    mode = 'auto'

    # Optional external MET CSV file (will be merged on the 'time' coordinate)
    # Set include_external_met_csv = False to disable.
    external_met_csv = r'D:\OneDrive - Universite de Liege\Documents\Science\Programmation\Projets Python\be-vie_tof_ec_analysis\L2_O3_files\BE-Lon_L2_2201-2605.csv'
    include_external_met_csv = False

    # Optional time window selection (applied to both HDF5 data and external MET)
    # Use day strings 'YYYY-MM-DD'. Keep None to disable a bound.
    start_time = '2022-04-08'
    end_time = '2026-04-12'

    if start_time or end_time:
        start_tag = start_time if start_time else ''
        end_tag = end_time if end_time else ''
        if start_tag and end_tag:
            window_tag = f"_{start_tag}_to_{end_tag}"
        elif start_tag:
            window_tag = f"_{start_tag}_to"
        else:
            window_tag = f"_to_{end_tag}"
        p = Path(output_nc)
        output_nc = str(p.with_name(f"{p.stem}{window_tag}{p.suffix}"))

    ds = concat_geddysoft_hdf5(
        results_folder=results_folder,
        prefix=prefix,
        output_nc=output_nc,
        mode=mode,
        start_time=start_time,
        end_time=end_time,
    )

    if 'time' in ds.coords and ds.sizes.get('time', 0) > 0:
        print(f"Time range (raw): {ds['time'].values[0]} -> {ds['time'].values[-1]}")

    ds = _apply_time_window(ds, start_time, end_time)
    ds = _drop_duplicate_times(ds)
    if 'time' in ds.coords and ds.sizes.get('time', 0) > 0:
        print(f"Time range (filtered): {ds['time'].values[0]} -> {ds['time'].values[-1]}")

    if include_external_met_csv and external_met_csv:
        met_ds = _read_external_met_csv(external_met_csv)
        met_ds = _apply_time_window(met_ds, start_time, end_time)
        met_ds = _drop_duplicate_times(met_ds)
        d1 = _count_duplicate_times(ds)
        d2 = _count_duplicate_times(met_ds)
        print(f"Duplicate time stamps (ds): {d1}; (external MET): {d2}")
        ds = xr.merge([ds, met_ds], join='left')

    output_path = Path(output_nc)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(output_path, engine='scipy', format='NETCDF3_64BIT')
    print(f"Wrote NetCDF: {output_path}")
