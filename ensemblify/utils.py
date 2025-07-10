import itertools

_SPACER = "§"


def equals_none(string:str):
    if string is None:
        return True
    string = string.strip().lower()
    if (string == "null") | (string == "none") | (string == "") | (string == "~"):
        return True
    return False


def get_option_variants(config:dict):
        """Returns a list with the cartesian product of all combinations of options, formatted as strings."""
        if config.get("params") is None:
            # Return empty string to avoid cmd error (necessary ?)
            return [""]
        
        variants = []
        base_options = []
        product_options = {}
        for key, val in config["params"].items(): 
            if equals_none(str(val)):
                # If val is none: simple flag
                base_options.append(key)
            elif isinstance(val, list):
                # If val is list: contains options
                product_options[key] = val
            else:
                # Val is assumed to be a single value
                base_options.append(key)
                base_options.append(str(val))

        product_keys = list(product_options.keys())
        for product_vals in itertools.product(*product_options.values()):
            variant = base_options.copy()
            for key, val in zip(product_keys, product_vals):
                variant.extend([key, str(val)])
            variants.append(variant)
        variants = [_SPACER.join(variant) for variant in variants]
        return variants