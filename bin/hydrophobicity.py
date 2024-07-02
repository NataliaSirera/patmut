def hydrophobicity(hydrophobicityDB, variants, absoluteFlag):
    if absoluteFlag:
        hydrophobicityARFF = dict(
            (
                variant,
                str(
                    abs(float(hydrophobicityDB[variant[0]]) - float(hydrophobicityDB[variant[-1]]))
                ),
            )
            for variant in variants
        )
    else:
        hydrophobicityARFF = dict(
            (
                variant,
                str(float(hydrophobicityDB[variant[0]]) - float(hydrophobicityDB[variant[-1]])),
            )
            for variant in variants
        )
    return hydrophobicityARFF
