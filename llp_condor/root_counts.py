"""Read every TTree's entry count, including trees in TFileService directories."""
import json
import sys


def tree_counts(filename):
    import ROOT
    ROOT.gROOT.SetBatch(True)
    f = ROOT.TFile.Open(str(filename), 'READ')
    if not f or f.IsZombie() or f.TestBit(ROOT.TFile.kRecovered):
        raise RuntimeError('ROOT file is missing, corrupt or recovered: ' + str(filename))
    counts = {}
    def walk(directory, prefix=''):
        # Consider the newest cycle of each object only.
        for name in sorted({key.GetName() for key in directory.GetListOfKeys()}):
            obj = directory.Get(name)
            path = prefix + name
            if obj.InheritsFrom('TDirectory'):
                walk(obj, path + '/')
            elif obj.InheritsFrom('TTree'):
                counts[path] = int(obj.GetEntries())
    try:
        walk(f)
    finally:
        f.Close()
    if not counts:
        raise RuntimeError('No TTrees found in ' + str(filename))
    return counts


if __name__ == '__main__':
    print(json.dumps(tree_counts(sys.argv[1]), sort_keys=True))
