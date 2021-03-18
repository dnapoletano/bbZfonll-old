class messages:
    def __init__(self):
        self.HEADER = '\033[95m'
        self.OKBLUE = '\033[94m'
        self.OKGREEN = '\033[92m'
        self.WARNING = '\033[93m'
        self.FAIL = '\033[91m'
        self.ENDC = '\033[0m'
        self.BOLD = "\033[1m"

    def disable(self):
        self.HEADER = ''
        self.OKBLUE = ''
        self.OKGREEN = ''
        self.WARNING = ''
        self.FAIL = ''
        self.ENDC = ''

    def print_intro(self):
        string='bbHFONLL: Calculation of the FONLL-XS for bb ->h'
        self.info(string)
        
    def print_outro(self):
        string=r'Please cite : Phys.Lett. B763 (2016) 190-196 and Phys.Lett. B751 (2015) 331-337, when using this code' 
        self.info(string)
        
    def header(self, msg):
        message = self.HEADER +msg +self.ENDC
        print(message)
        
    def infog(self, msg):
        message = self.OKGREEN + msg + self.ENDC
        print(message)

    def info(self, msg):
        message = self.OKBLUE + msg + self.ENDC 
        print(message)

    def warn(self, msg):
        message = self.WARNING + msg + self.ENDC
        print(message)

    def err(self, msg):
        message = self.FAIL + msg + self.ENDC
        print(self.FAIL + msg + self.ENDC)
