def substitutionMatrix(substitutionMatrixDB, variants):
    substitutionMatrixARFF = dict(
        (
            variant,
            str(substitutionMatrixDB[(variant[0], variant[-1])])
            if (variant[0], variant[-1]) in substitutionMatrixDB
            else str(substitutionMatrixDB[(variant[-1], variant[0])]),
        )
        for variant in variants
    )
    return substitutionMatrixARFF
