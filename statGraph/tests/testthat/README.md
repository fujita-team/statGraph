# Contributing to the codebase
## Guidelines

All test files should be named "test-{*method_name*}.R". 

In order to test a block of code, use the function RunTest, for example:

    RunTest( {
        *Your Code Here*
    })

Alternative, set your code as a variable, and invoke RunTest by calling this variable, ex:

    test_code <- { 
        *Your Code Here* 
    }
    RunTest(test_code)`

Do **NOT** run code outside of the method RunTest, otherwise consistency of test results cannot be assured.

For indications on how to verify certain conditions, please read: https://testthat.r-lib.org/reference/

Do **NOT** remove the file *helper-methods.R*
