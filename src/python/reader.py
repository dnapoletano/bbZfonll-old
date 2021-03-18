import yaml
import numpy as np

class reader:
    def __init__(self,namefile):
        var_names, var_vals = [], []
        with open(namefile,'r') as stream:
            try:
                docs=yaml.load_all(stream)
                for doc in docs:
                    for k,v in doc.items():
                        print(k,"-->",v)
                        var_names.append(k)
                        var_vals.append(v)
            except yaml.YAMLError as exc:
                    print(exc)
        self.var_read=np.array([var_names,var_vals])
        self.variables=np.zeros(len(var_names),
                                dtype = {'names': ['Var Name', 'Var Value'],
                                     'formats': ['U20', 'U20']} )
        self.variables['Var Name'] = self.var_read[0]
        self.variables['Var Value' ] = self.var_read[1]
                
    def get_variables(self):
        return self.variables
    
    def get_variable(self,variable,default):
        if(np.any(self.variables['Var Name']==str(variable))):
            index=list(self.variables['Var Name']).index(str(variable))
            return self.variables['Var Value'][index]
        else:
            return default
    
    def get_d_variable(self,variable,default):
        if(np.any(self.variables['Var Name']==str(variable))):
            index=list(self.variables['Var Name']).index(str(variable))
            return float(self.variables['Var Value'][index])
        else:
            return default
    
    def get_i_variable(self,variable,default):
        if(np.any(self.variables['Var Name']==str(variable))):
            index=list(self.variables['Var Name']).index(str(variable))
            return int(self.variables['Var Value'][index])
        else:
            return default
    
    def get_str_variable(self,variable,default):
        if(np.any(self.variables['Var Name']==str(variable))):
            return self.get_variable(variable,default)
        else:
            return default
