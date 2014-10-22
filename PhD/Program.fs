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
    elitism = 0.1; mutate = 0.8; spawn = 0.0; selection = 0.4; newBlood = 0.05;
    creator = testCreate}
     
let evolve ga popSize steps = 
    let r = System.Random()
    Seq.fold (fun (s:'a[]) _ -> printfn "%f" (ga.fitness s.[0]); runGA r ga s) 
        (Array.init popSize (fun _ -> ga.creator r)) 
        (seq {1..steps})

[<EntryPoint>]
let main argv = 
    let m = {
        stem=stemScore; internalDangle = internalDangleScore; extenalDangle = extenalDangleScore; lbulge = bulgeScore;
        internalLoop = internalLoopScore; hairpin=hairpinLoopScore; rbulge = bulgeScore; externalStemBonus = externalStemBonus; interalStemBonus = interalStemBonus} 
    let rna = System.Console.ReadLine() |> RNAPrimary.parse |> Seq.toArray
    let model = BasicModel.makeModel rna (BasicModel.Paramaters.Random())
    printfn "%A" (MFE.zuker rna model)
    System.Console.Read()
    0 // return an integer exit code
