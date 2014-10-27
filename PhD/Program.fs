open RNASecondary
open MFE
open GA

let stemScore i j k = 1.0
let internalLoopScore i k j l = 1.0
let bulgeScore i j = 1.0
let hairpinLoopScore i j = 1.0
let internalDangleScore i j = 1.0
let extenalDangleScore i j= 1.0
let externalStemBonus i j = 1.0
let interalStemBonus i j = 1.0

let testBreeder _ a b = (a+b)/2.0
let testMutate (r:System.Random) i = i + (r.NextDouble() - 0.5)
let testFitness i = System.Math.Abs(12.32 - i)
let testCreate (r:System.Random) = r.NextDouble()
let testGA = {breeder = testBreeder; mutator = testMutate; fitness = testFitness;
    elitism = 0.1; mutate = 0.8; selection = 0.4; newBlood = 0.05;
    creator = testCreate}
     
let evolve ga popSize steps = 
    let r = System.Random()
    Seq.fold (fun (s:'a[]) _ -> printfn "%f" (ga.fitness s.[0]); runGA r ga s) 
        (Array.init popSize (fun _ -> ga.creator r)) 
        (seq {1..steps})

let r = System.Random()
let randomSequence sz = 
    let b = "augc"
    new string [| for i in 1..sz -> b.[r.Next(3)]|]
let testSequence s = 
    let rna = s |> RNAPrimary.parse |> Seq.toArray
    let par = FullModel.Parameters.Random()
    let model = FullModel.makeModel rna par
    let ssScore = RNASecondary.parseSurfaces >> RNASecondary.scoreExternalLoop model
    let best = rna |> RNASecondary.validSecondaryStructures |> Seq.maxBy ssScore
    let bests = ssScore best
    let zuks = MFE.zuker rna model
    if System.Math.Abs(bests-zuks) > 0.00001 then
        printfn "%A vs %A \n\n %A \n\n %A \n\n %A" bests zuks best rna par
    else
        ()

open BasicModel
let trainBasicModel rnas sstructs = 
    let testRNAs = Seq.zip (Seq.map (RNAPrimary.parse>>Seq.toArray) rnas) (Seq.map RNASecondary.parseSurfaces sstructs)
    let ga = 
        BasicModel.GATrainer
            {
                mutationRate = 0.3;
                mutationRange = 1.0;
                crossoverRate = 0.5;
                testRNAs = testRNAs;
            }
    let r = System.Random()
    Seq.fold 
        (fun (s:Genome[]) t -> 
            printfn "Generation #%d: \n %A \n\n" t s.[0]
            GA.runGA r ga s
        )
        (Array.init 100 (fun _ -> ga.creator r))
        (seq {1..250})

let fullModelGA testRNAs =
    let fitness = FullModel.calcFitness testRNAs 
    {
        GA.breeder = FullModel.breed fitness 0.4;
        GA.mutator = FullModel.mutate fitness 0.5 1.0;
        GA.creator = 
            fun _ -> 
                let p = FullModel.Parameters.Random()
                { FullModel.parameters = p; FullModel.fitness = fitness p };
        GA.fitness = fun a -> a.fitness;
        elitism = 0.05;
        newBlood = 0.05;
        mutate = 0.4;
        selection = 0.45;
    }

let trainFullModel rnas sstructs = 
    let testRNAs = 
        Seq.zip 
            (Seq.map (RNAPrimary.parse>>Seq.toArray) rnas) 
            (Seq.map RNASecondary.parseSurfaces sstructs)

    let r = System.Random()
    let ga = fullModelGA testRNAs
    Seq.fold 
        (fun (s:FullModel.Genome[]) t -> 
            printfn "Generation #%d: \n %A \n\n" t s.[0]
            GA.runGA r ga s
        )
        (Array.init 100 (fun _ -> ga.creator r))
        (seq {1..250})



[<EntryPoint>]
let main argv = 
    for i in 1..10000 do
        testSequence (randomSequence 3)
    System.Console.Read() |> ignore
    0 // return an integer exit code
