module DP

let memoize baseCases (fn:'a->'b) = 
    let d = System.Collections.Generic.Dictionary<'a, 'b>()
    for (bi, bo) in baseCases do
        d.Add(bi, bo)
    (fun inp -> 
        if d.ContainsKey(inp) then d.[inp] 
        else 
            let v = fn inp
            d.Add(inp, v)
            v)
