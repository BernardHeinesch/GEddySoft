import math


def remove_entries_based_on_indices(data, nan_indices, exceptions=[]):
    if isinstance(data, dict):
        for key, value in data.items():
            if key in exceptions:
                continue
            if isinstance(value, list):
                data[key] = [v for i, v in enumerate(value) if i not in nan_indices]
            elif isinstance(value, dict):
                remove_entries_based_on_indices(value, nan_indices)
    elif isinstance(data, list):
        return [v for i, v in enumerate(data) if i not in nan_indices]


def clean_results(results):
    """
        Function to remove from results all unuseful lines, i.e. all lines for which no eligible sonic data where available
        (e.g. incomplete sonic files). The identification of an unuseful lines is made on Nan values in the wsh list in MET.

        parameters
        ----------
        results: dict, final dictionnary of GEddySoft

        returns
        -------
        results: dict, final dictionnary of GEddySoft, with unuseful entries eliminated

        comments
        --------
        Written by B. Heinesch.
        University of Liege, Gembloux Agro-Bio Tech.
    """

    met = results.get('MET', {})

    # Identify indices with NaN values in 'wsh'
    nan_indices = [i for i, value in enumerate(met['wsh']) if math.isnan(value)]

    # Remove corresponding entries from 'time'
    results['time'] = [time for i, time in enumerate(results['time']) if i not in nan_indices]

    # Remove entries from MET based on nan_indices
    remove_entries_based_on_indices(met, nan_indices, exceptions=['qaqc'])

    # Remove entries from MET.qaqc based on nan_indices
    remove_entries_based_on_indices(met['qaqc'], nan_indices)

    # Remove entries from each group in IRGA based on nan_indices
    tracer = results.get('IRGA', {})
    for group_key, group_value in tracer.items():
        if isinstance(group_value, dict):
            remove_entries_based_on_indices(group_value, nan_indices, exceptions=['name', 'qaqc'])
            if 'qaqc' in group_value:
                remove_entries_based_on_indices(group_value['qaqc'], nan_indices)
        elif isinstance(group_value, list):
            tracer[group_key] = remove_entries_based_on_indices(group_value, nan_indices)

    # Remove entries from each group in TRACER based on nan_indices
    tracer = results.get('TRACER', {})
    for group_key, group_value in tracer.items():
        if isinstance(group_value, dict):
            remove_entries_based_on_indices(group_value, nan_indices, exceptions=['name', 'qaqc'])
            if 'qaqc' in group_value:
                remove_entries_based_on_indices(group_value['qaqc'], nan_indices)
        elif isinstance(group_value, list):
            tracer[group_key] = remove_entries_based_on_indices(group_value, nan_indices)

    return (results)
