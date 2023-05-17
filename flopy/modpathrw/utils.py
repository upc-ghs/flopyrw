'''
General utilities 
'''


def count_instances(method):
    '''
    Decorator for counting class instances. 
    Rely on the class having a COUNTER variable
    '''

    def wrapper(self, *args, **kw):
        self.__class__.COUNTER += 1
        return method(self, *args, **kw)

    return wrapper
