# DO NOT DELETE THIS FILE!
# THIS FILE SETS UP THE RANDOM SEED FOR THE REMAINING TESTS!
#
# DON'T RUN MANUALLY!
# THIS SCRIPT IS INVOKED AUTOMATICALLY BY devtools::test()
# Run 'Rscript tests.R' to perform tests or 'Rscript build.r' to build the project.
# For Documentation on how to create tests, please refer to: https://cran.r-project.org/web/packages/testthat/testthat.pdf


generateSeed <- function(string){
    prime <- 2428095424619
    result <- 0
    for(char in strsplit(string, "")[[1]]){
        result <- (result + utf8ToInt(char)) %% prime
    }
    return(result)
}


poem <-
    "
                Two-Headed Calf

    Tomorrow when the farm boys find this
    freak of nature, they will wrap his body
    in newspaper and carry him to the museum.

    But tonight he is alive and in the north
    field with his mother. It is a perfect
    summer evening: the moon rising over
    the orchard, the wind in the grass. And
    as he stares into the sky, there are
    twice as many stars as usual."


RunTest <- function(Name, Code, seed = generateSeed(poem)){
    set.seed(seed)
    if(is.null(VALID_TESTS)){
        eval.parent(Code)
        return()
    }
    if(VALID_TESTS$mode == "ALL")
    {
        eval.parent(Code)
        return()
    }
    if(VALID_TESTS$mode == "WHITELIST")
    {
        if(Name %in% VALID_TESTS$list)
        {
            eval.parent(Code)
            return()
        }
    }
    if(VALID_TESTS$mode == "BLACKLIST")
    {
        if(! Name %in% VALID_TESTS$list)
        {
            eval.parent(Code)
            return()
        }
    }
}

