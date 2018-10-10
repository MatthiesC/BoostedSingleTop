channels = ["Electron", "Muon"]
ptbins   = ["200to300", "300to400", "400to600", "600to1200"]
nbtags   = ["1b", "2b"]

steerTemplate = "ThetaTemplate.steer"

for stage in ["Prefit", "Postfit"]:
    for c in channels:
        for p in ptbins:
            for nb in nbtags:
                steerNew = "Theta"+c+nb+p+stage+".steer"
                with open(steerTemplate, "rt") as fin:
                    with open(steerNew, "wt") as fout:
                        for line in fin:
                            line_new = line.replace("<inputname>",("input" if stage == "Prefit" else "output")+"_"+c+"_"+nb+"_"+p)
                            line_new = line_new.replace("<outputname>",c+stage)
                            fout.write(line_new)
