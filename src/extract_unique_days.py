from datetime import datetime

def extract_unique_days(filelist, prefix=""):
    """
    Extract unique days from filenames, keeping only year, month, day
    in the same format as in the filenames.

    Parameters
    ----------
    filelist: list of str
        List of filenames
    prefix: str
        Optional prefix to strip before parsing

    Returns
    -------
    sorted list of unique day strings, in same Y/M/D format as in filenames
    """

    days = set()

    # known formats: (format string, length of date part, length of Y/M/D part)
    formats = [
        ("%Y-%m-%dT%H%M%S", 17, 10),        # 2023-07-06T000000 → Y-M-D length=10
        ("%Y_%m_%d__%H_%M_%S", 20, 10),     # 2023_07_06__00_00_00 → Y_M_D length=10
        ("%Y%m%d%H%M", 12, 8),              # 202411010000 → YMD length=8
    ]

    for filename in filelist:
        name = filename

        # remove prefix if present
        if prefix and name.startswith(prefix):
            name = name[len(prefix):]

        # remove extension
        name = name.split('.')[0]

        parsed = False
        for fmt, total_length, ymd_length in formats:
            try:
                # validate datetime
                _ = datetime.strptime(name[:total_length], fmt)
                # store only the Y/M/D substring
                days.add(name[:ymd_length])
                parsed = True
                break
            except ValueError:
                pass

        if not parsed:
            raise ValueError(f"Unrecognized filename format: {filename}")

    return sorted(days)
