class IncrementingDictionary(dict):
    """
    A dictionary like object that adds new keys to the dictionary with an incremented index as its value.
    """

    def __init__(self, *args, **kwargs):
        super(IncrementingDictionary, self).__init__(*args, **kwargs)

    def __setitem__(self, key, value):
        if key in self:
            raise KeyError('The value associated with a key in a ConstraintMap can only be set once')
        else:
            super(IncrementingDictionary, self).__setitem__(key, len(self))

    def __getitem__(self, key):
        try:
            return super(IncrementingDictionary, self).__getitem__(key)
        except KeyError:
            self.__setitem__(key, None)
            return super(IncrementingDictionary, self).__getitem__(key)