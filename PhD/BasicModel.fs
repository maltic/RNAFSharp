module BasicModel

open RNASecondary

type Paramaters = 
    {
        stacking : float[,];
        unpaiedExteral : float[];
        unpairedInternal : float[];
        lbulge : float[];
        rbulge : float[];
        hairpin : float[];
        internalLoop : float[];
        externalStemM : float;
        externalStemC : float;
        internalStemM : float;
        internalStemC : float;
    }
    static member Empty = 
        let b = RNAPrimary.BASES
        {
            stacking = Array2D.zeroCreate b b;
            unpaiedExteral = Array.zeroCreate b;
            unpairedInternal = Array.zeroCreate b;
            lbulge = Array.zeroCreate b;
            rbulge = Array.zeroCreate b;
            hairpin = Array.zeroCreate b;
            internalLoop = Array.zeroCreate b;
            externalStemM = 0.0;
            externalStemC = 0.0;
            internalStemM = 0.0;
            internalStemC = 0.0;
        }
    static member Random() = 
        let r = System.Random()
        let rand _ = r.NextDouble()
        let b = RNAPrimary.BASES
        {
            stacking = Array2D.init b b (fun _ _ -> r.NextDouble());
            unpaiedExteral = Array.init b rand;
            unpairedInternal = Array.init b rand;
            lbulge = Array.init b rand;
            rbulge = Array.init b rand;
            hairpin = Array.init b rand;
            internalLoop = Array.init b rand;
            externalStemM = rand();
            externalStemC = rand();
            internalStemM = rand();
            internalStemC = rand();
        }

let stackingScore (param:float[,]) (rna:RNAPrimary.Base[]) i j k =
    Seq.sum (seq { for d in 0..k-1 -> param.[int rna.[i+d], int rna.[j-d]] })
let loopScore (param:float[]) (rna:RNAPrimary.Base[]) i j = 
    Seq.sum (seq { for k in i..j -> param.[int rna.[k]] })

let makeModel (rna:RNAPrimary.Base[])  paramaters = 
    {
        stem = stackingScore paramaters.stacking rna;
        internalLoop = 
            let f = loopScore paramaters.internalLoop rna
            fun i j k l -> f i j + f k l;
        lbulge = loopScore paramaters.lbulge rna;
        rbulge = loopScore paramaters.rbulge rna;
        hairpin = loopScore paramaters.hairpin rna;
        internalDangle = loopScore paramaters.unpairedInternal rna;
        extenalDangle = loopScore paramaters.unpaiedExteral rna;
        externalStemBonus = 
            fun i j -> paramaters.externalStemM * float (j-i+1) + paramaters.externalStemC;
        interalStemBonus = 
            fun i j -> paramaters.internalStemM * float (j-i+1) + paramaters.internalStemC;
    }

type Genome = 
    {
        fitness : double;
        parameters : Paramaters;
    }

type GASettings = 
    {
        mutationRate : float;
        mutationRange : float;
        crossoverRate : float;
        computeFitness : Paramaters -> double;
    }

let calcFitness testRNAs parameters =
    let helper (rna,ss) = 
        let genomeModel = makeModel rna parameters
        let mfe = MFE.zuker rna genomeModel
        let actual = RNASecondary.scoreExternalLoop genomeModel ss
        mfe - actual
    Seq.sumBy helper testRNAs

let mutate settings (r:System.Random) g = 
    let maybeMutate a = 
        if r.NextDouble() > settings.mutationRate then
            a
        else
            a + ((r.NextDouble() * settings.mutationRange) - settings.mutationRange/2.0)
    let gparams = g.parameters
    let p = 
        {
            stacking = Array2D.map maybeMutate gparams.stacking;
            unpaiedExteral = Array.map maybeMutate gparams.unpaiedExteral;
            unpairedInternal = Array.map maybeMutate gparams.unpairedInternal;
            lbulge = Array.map maybeMutate gparams.lbulge;
            rbulge = Array.map maybeMutate gparams.rbulge;
            hairpin = Array.map maybeMutate gparams.hairpin;
            internalLoop = Array.map maybeMutate gparams.internalLoop;
            externalStemM = maybeMutate gparams.externalStemM;
            externalStemC = maybeMutate gparams.externalStemC;
            internalStemM = maybeMutate gparams.internalStemM;
            internalStemC = maybeMutate gparams.internalStemC;
        }
    {
        fitness = settings.computeFitness p;
        parameters = p;
    }
