def freqGaps(msaDBT, alignedQuery, variants):
    freqGaps = [
        col.count("-") / float(len(msaDBT[0]))
        for pos, col in enumerate(msaDBT)
        if alignedQuery[pos] != "-"
    ]
    freqGapsARFF = dict((variant, freqGaps[int(variant[1:-1]) - 1]) for variant in variants)
    return freqGapsARFF
