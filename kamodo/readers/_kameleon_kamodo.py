from kamodo import Kamodo, kamodofy
from . import kameleon_gateway
import numpy as np

import execnet

class Kameleon(Kamodo):
    def __init__(self, filename, python_path, *variables, **kwargs):
        
        self._current_file = filename
        self._python_path = python_path
        self.make_gateway()
        
        self.open_file(filename, *variables)
        
        super(Kameleon, self).__init__(**kwargs)
        
        self.make_interpolators()
        
        
    def get_coords(self):
        """Extracts coordinate list from grid_system_1 metadata:
            STRING: grid_system_1: [X,Y,Z]
        '"""
        coords = self._metadata['grid_system_1'].split(':')[-1].strip()
        coords = ''.join(c for c in coords if c not in '[]') # remove brackets
        coords = [c.strip().lower() for c in coords.split(',')]
        return coords
        
    def make_gateway(self):
        """Create a gateway into python installed with kameleon"""
        gateway_spec = 'popen//dont_write_bytecode//python={}'.format(self._python_path)
        self._gw = execnet.makegateway(gateway_spec)
        self._ch = self._gw.remote_exec(kameleon_gateway)
        
    def open_file(self, file_name, *variables):
        args = [file_name] + list(variables)
        cmd_args = ", ".join(["'{}'".format(s) for s in args])
        cmd = "initialize({})".format(cmd_args)
        try:
            self._ch.send(cmd) # execute func-call remotely
        except self._ch.RemoteError as m:
            print(m)
            
        try:
            self._metadata = self._ch.receive()
            for variable_name in variables:
                assert variable_name in self._metadata['variables']
            # sort variable metadata in the same order as input
            self._metadata['variables'] = {variable : self._metadata['variables'][variable] \
                                           for variable in variables}
            
            self._current_file = file_name
        except:
            print('could not open file')
            print('offending command:\n\t{}'.format(cmd))
            print(self._metadata)
        print('{} opened'.format(file_name))
        
        
    def interpolate(self, variable_name, c0, c1, c2):

        self._ch.send("interpolate('{}',{},{},{})".format(variable_name, c0, c1, c2))

        return self._ch.receive()

    def make_interpolators(self):
        for variable, metadata in self._metadata['variables'].items():
            coords = self.get_coords()
            units = metadata['units']
            exec_str = """def interpolate({c0}vec):
                points = {c0}vec
            
                if type(points) == np.ndarray:
                    points = points.tolist()
                package = dict(
                    variable = '{varname}', 
                    points = points)
                self._ch.send(package)
                results = self._ch.receive()
                if type({c0}vec) == np.ndarray:
                    return np.array(results)
                else:
                    return results
            """.format(varname = variable,
                      c0 = coords[0],
                      c1 = coords[1],
                      c2 = coords[2])
            d = {'self':self, 'np':np}
            exec(exec_str, d)
            interp_func = d['interpolate']
            self[variable] = kamodofy(interp_func, units = units)