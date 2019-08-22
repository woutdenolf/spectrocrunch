# -*- coding: utf-8 -*-

class defaultdict(dict):
    """Dictionary with a lambda function on the key value as default value
    """

    def __init__(self, *arg, **kw):
        super(defaultdict, self).__init__(*arg, **kw)
        self.lambdafunc = lambda key: key

    def setdefaultfactory(self, lambdafunc):
        self.lambdafunc = lambdafunc

    def __getitem__(self, key):
        if key in self:
            return super(defaultdict, self).__getitem__(key)
        else:
            return self.lambdafunc(key)
