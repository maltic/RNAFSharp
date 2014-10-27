module FullModel


type Parameters = 
    {
        stem : Stem.Parameters;
        HP: Hairpin.Parameters;
        bulge : Bulge.Parameters;
        IL : InternalLoop.Parameters;
        ML : MultiLoop.Parameters;
    }
    static member Random() = 
        {
            stem = Stem.Parameters.Random ();
            HP = Hairpin.Parameters.Random ();
            bulge = Bulge.Parameters.Random ();
            IL = InternalLoop.Parameters.Random ();
            ML = MultiLoop.Parameters.Random ();
        }

let makeModel rna p = 
    {
        RNASecondary.stem = Stem.score p.stem rna;
        RNASecondary.internalLoop = InternalLoop.score p.IL rna;
        RNASecondary.lbulge = Bulge.score p.bulge rna;
        RNASecondary.rbulge = Bulge.score p.bulge rna;
        RNASecondary.hairpin = Hairpin.score p.HP rna;
        RNASecondary.internalDangle = MultiLoop.scoreInternalDangle p.ML rna;
        RNASecondary.interalStemBonus = MultiLoop.scoreInternalStem p.ML rna;
        RNASecondary.externalDangle = MultiLoop.scoreExternalDangle p.ML rna;
        RNASecondary.externalStemBonus = MultiLoop.scoreExternalStem p.ML rna;
    }

type Genome =
    {
        parameters : Parameters;
        fitness : float;
    }

let calcFitness testRNAs parameters =
    let helper (rna,ss) = 
        let genomeModel = makeModel rna parameters
        let mfe = MFE.zuker rna genomeModel
        let actual = RNASecondary.scoreExternalLoop genomeModel ss
        if System.Double.IsNaN(mfe) then printfn "mfe is NaN"
        elif System.Double.IsInfinity(mfe) then printfn "mfe is inf"
        elif System.Double.IsNaN(actual) then printfn "actual is NaN"
        elif System.Double.IsInfinity(actual) then printfn "actual is inf"
        mfe - actual
    Seq.sumBy helper testRNAs

let mutate fitness rate range (r:System.Random) genome = 
    let newParams = 
        {
            stem = Stem.mutate rate range r genome.parameters.stem;
            HP = Hairpin.mutate rate range r genome.parameters.HP;
            bulge = Bulge.mutate rate range r genome.parameters.bulge;
            IL = InternalLoop.mutate rate range r genome.parameters.IL;
            ML = MultiLoop.mutate rate range r genome.parameters.ML;
        }
    {
        parameters = newParams;
        fitness = fitness newParams;
    }

let breed fitness rate (r:System.Random) a b = 
    let newParams = 
        {
            stem = Stem.breed rate r a.parameters.stem b.parameters.stem;
            HP = Hairpin.breed rate r a.parameters.HP b.parameters.HP;
            bulge = Bulge.breed rate r a.parameters.bulge b.parameters.bulge;
            IL = InternalLoop.breed rate r a.parameters.IL b.parameters.IL;
            ML = MultiLoop.breed rate r a.parameters.ML b.parameters.ML;
        }
    {
        parameters = newParams;
        fitness = fitness newParams;
    }

