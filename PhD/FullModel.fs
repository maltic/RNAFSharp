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

let makeModel p rna = 
    {
        RNASecondary.stem = Stem.score p.stem rna;
        RNASecondary.internalLoop = InternalLoop.score p.IL rna;
        RNASecondary.lbulge = Bulge.score p.bulge rna;
        RNASecondary.rbulge = Bulge.score p.bulge rna;
        RNASecondary.hairpin = Hairpin.score p.HP rna;
        RNASecondary.internalDangle = MultiLoop.scoreInternalDangle p.ML rna;
        RNASecondary.interalStemBonus = MultiLoop.scoreInternalStem p.ML rna;
        RNASecondary.extenalDangle = MultiLoop.scoreExternalDangle p.ML rna;
        RNASecondary.externalStemBonus = MultiLoop.scoreExternalStem p.ML rna;
    }

