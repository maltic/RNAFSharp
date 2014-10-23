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
    let params = BasicModel.Paramaters.Random()
    let model = BasicModel.makeModel rna params
    let ssScore = RNASecondary.parseSurfaces >> RNASecondary.scoreExternalLoop model
    let best = rna |> RNASecondary.validSecondaryStructures |> Seq.maxBy ssScore
    let bests = ssScore best
    let zuks = MFE.zuker rna model
    if System.Math.Abs(bests-zuks) > 0.00001 then
        printfn "%A vs %A \n\n %A \n\n %A \n\n %A" bests zuks best rna params
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
            printfn "Generation #%d: %f \n\n %A \n\n" t (ga.fitness s.[0]) s.[0]
            GA.runGA r ga s
        )
        (Array.init 25 (fun _ -> ga.creator r))
        (seq {1..100})


[<EntryPoint>]
let main argv = 
    trainBasicModel ["GAUCACUUCGGUGAUCACUUCGGU"; "GGCCAGAUUGAGCCUGGGAGCUCUCUGGCC"] ["....((((((((....))))))))"; "(((((((..((((......)))))))))))"]
    System.Console.Read() |> ignore
    0 // return an integer exit code
