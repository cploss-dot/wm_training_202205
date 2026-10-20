#!/usr/bin/env nextflow

/* 
this is now deprecated in nextflow since DSL2 is default
*/
nextflow.enable.dsl = 2

/*
These are just basic channel manipulations. But one thing to consider
is they don't execute in the order they appear below. Sometime .count()
is first, sometimes its last, sometimes its in the middle somewhere.

Learn more about Operators at: https://www.nextflow.io/docs/edge/operator.html#
*/

/*
added a workflow block -> otherwise causes an error
*/
workflow {
    Channel
        .from( 1, 2, 3, 4, 5 )
        .view{ num -> "$num" }

    Channel
        .from( 1, 2, 3, 4, 5 )
        .count()
        .view{ num -> "Output of .count(): $num" }

    Channel
        .from( 1, 2, 3, 4, 5 )
        .map { it -> [it, it*it] }
        .view { num, sqr -> "Square of: $num is $sqr" }

    Channel
        .from( 1, 2, 3, 4, 5 )
        .toList()
        .view { num -> ".toList(): $num" }

    Channel
        .from( 1, 2, 3, 4, 5 )
        .sum( { it } )
        .view { ".sum(): $it" }

    Channel
        .from( 1, 2, 3, 4, 5 )
        .min( { it } )
        .view { ".min(): $it" }

    Channel
        .from( 1, 2, 3, 4, 5 )
        .randomSample( 2 )
        .view{ ".randomSample(): $it" }

    // added max operator example
    Channel
        .of( 8, 6, 10, 5, -9 )
        .max()
        .view { v -> "Max value is $v" }

    // added collect opetator example
    Channel
        .of( -1, 2, -3, 4 )
        .collect()
        .view()
}