﻿module Hairpin

type Parameters =
    {
        a : float;
        b : float;
        c : float;
        values : float[];
    }
    static member Random() = 
        let r = System.Random()
        let rand _ = r.NextDouble()
        let b = RNAPrimary.BASES
        {
            a = rand();
            b = rand();
            c = rand();
            values = Array.init b rand
        }

let score p (rna:RNAPrimary.Base[]) i j = 
    ModelFunctions.rnaRange rna i j 
        |> ModelFunctions.jacobsonStockmayer p.a p.b p.c p.values