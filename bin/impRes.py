def impRes(impResDB, variants):
    impResARFF = dict(
        (variant, "True" if variant[1:-1] in impResDB else "False") for variant in variants
    )
    return impResARFF
