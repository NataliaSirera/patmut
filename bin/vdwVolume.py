def vdwVolume(vdwVolumeDB, variants, absoluteFlag):
    if absoluteFlag:
        vdwVolumeARFF = dict(
            (variant, str(abs(float(vdwVolumeDB[variant[0]]) - float(vdwVolumeDB[variant[-1]]))))
            for variant in variants
        )
    else:
        vdwVolumeARFF = dict(
            (variant, str(float(vdwVolumeDB[variant[0]]) - float(vdwVolumeDB[variant[-1]])))
            for variant in variants
        )
    return vdwVolumeARFF
