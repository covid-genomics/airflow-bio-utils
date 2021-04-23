def resolve_callable(val, *args, **kwargs):
    if callable(val):
        return val(*args, **kwargs)
    return val
