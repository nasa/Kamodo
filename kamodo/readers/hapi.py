import urllib, json
from kamodo import Kamodo, kamodofy
from collections import OrderedDict

def kamodofy_column(df, c):
    @kamodofy(units = c.split('[')[-1].split(']')[0])
    def f(t = df.index):
        return df[c].reindex(t).interpolate('linear')
    return f


def get_info(hapi_url, hapi_id):
    query = '{}/info?id={}'.format(hapi_url, hapi_id)
    response = urllib.urlopen(query)
    data = json.loads(response.read())
    return data

def get_date_range(hapi_url, hapi_id):
    info = get_info(hapi_url, hapi_id)
    start_date = info['startDate']
    end_date = info['stopDate']
    return start_date, end_date

def get_data_url(hapi_url, hapi_id, start_date, end_date):
    return hapi_url + '/data?id={}&time.min={}&time.max={}&format=json'.format(hapi_id, start_date, end_date)

def load_data(data_url, fill_values):
    response = urllib.urlopen(data_url)
    data = json.loads(response.read())
    names = [d['name'] for d in data['parameters']]
    units = [d['units'] for d in data['parameters']]
    variables = ['{}[{}]'.format(*v) for v in zip(names, units)]
    # before you set the index, apply fill values
    df = pd.DataFrame(data['data'])
    df.columns = fill_values.keys()
    for c in df.columns:
        try:
            df[c][df[c] == fill_values[c]] = np.NaN
        except:
            print('skipped filling of {}'.format(c))
#     df = pd.DataFrame(data['data'], columns = variables).set_index(variables[0])
    df.columns = variables
    df.set_index(variables[0], inplace = True)
    df.index = pd.to_datetime(df.index)
    return df

def get_catalog(hapi_url):
    response = urllib.urlopen(hapi_url + '/catalog')
    data = json.loads(response.read())
    catalog = pd.DataFrame(data['catalog'])
    return catalog

def get_fill_values(data):
    fill_vals = []
    for param in data['parameters']:
        try:
            fill_vals.append((param['name'], pd.to_numeric(param['fill'])))
        except:
            fill_vals.append((param['name'], np.NaN))
    return OrderedDict(fill_vals)

class Hapi(Kamodo):
    def __init__(self, hapi_url, hapi_id, start_date = None, stop_date = None, *args, **kwargs):
        super(Hapi, self).__init__(*args, **kwargs)
        
        start, end = get_date_range(hapi_url, hapi_id)
        if start_date is None:
            start_date = start
            
        if stop_date is None:
            stop_date = end
            
        data_url = get_data_url(hapi_url, hapi_id, start_date, stop_date)
        fill_values = get_fill_values(get_info(hapi_url, hapi_id))
        
        data = load_data(data_url, fill_values)
        
        for c in data.columns:
            self[c] = kamodofy_column(data, c)