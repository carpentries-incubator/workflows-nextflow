---
layout: reference
---

## Groovy syntax

Here are some fundamental
concepts of the Groovy language.

- Single line comments (lines not meant to be executed) start
with `//`. Multi-line comments are nested between `/*` and `*/` tags.
    ~~~
    // This is a single line comment. Everything after the // is ignored.
    /*
        Comments can also
        span multiple
        lines.
    */
    ~~~
    {: .language-groovy}
- Variables are assigned using `=` and can have any value. Variables used
inside a double quoted string are prefixed with a `$` to denote the
variable should be interpolated (replace the variable with the variable's
value). A variable is otherwise just used by it's
name.
    ~~~
    myvar = 1                                     // Integer
    myvar = -3.1499392                            // Floating point number
    myvar = false                                 // Boolean
    myvar = "Hello world"                         // String
    myvar = new java.util.Date()                  // Object - Abstract data structure
    message = "The file $file cannot be found!"   // Variable used inside a string.
    println message                               // Print variable
    ~~~
    {: .language-groovy}
- Lists (also known as arrays) are defined using the square bracket `[]` notation.
    ~~~
    emptylist = []                      // an empty list
    // Lists can contain duplicates, and the values can be of any type.
    mixedList = [1, 1, 'string', true, null, 5 as byte]
    mixedList.add("new value")          // adds "new value" to the end of mixedList
    println mixedList.size()            // prints the size
    println mixedList[0]                // prints 1
    ~~~
    {: .language-groovy}
- Maps (also known as associative arrays) are defined using the `[:]` literal. They associate a unique string with a value, and are commonly referred to as key-value pairs.
    ~~~
    emptyMap = [:]                      // an empty map
    mymap = [ name : "Steve", age: 43, likes: ['walks','cooking','coding']]
    mymap.age = 44                      // values can be accessed using map property notation.
    println mymap['name']               // alternatively values can be accessed using quoted keys.
    ~~~
    {: .language-groovy}
- Groovy supports common control structures such as if/else tests,
switch/case tests, for loops, and while loops.
    ~~~
    // if / else test
    x = false
    if ( !x ) {
        return 1
    } else {
        return 0
    }

    // switch / case test
    switch (x) {
        case "found foo":
            result = "Got foo"
            // continues to following cases

        case Number:
            result = "Number"
            break          // stops at this case

        default:
            result = "default"
    }

    // for loop (other variants exist)
    String message = ''
    for (int i = 0; i < 5; i++) {
        message += 'Hi '
    }

    // while loop
    y = 5
    while ( y-- > 0 ){
        println "Are we there yet?"
    }
    println "We've arrived!"
    ~~~
    {: .language-groovy}
- Closures are open, anonymous, blocks of code that can take arguments,
return a value and be assigned to a variable. A closure definition
follows the syntax `{ [closureParameters ->] statements }`. A closure
is evaluated when it is called, returning the result of the last statement
in the closure.
    ~~~
    // Find all elements > 1
    above_one = [1, 2, 3].findAll { it -> it > 1 }

    // Closure in a string
    message = "The file ${ file.getName() } cannot be found!"
    ~~~
    {: .language-groovy}

Groovy is very syntax-rich and supports many more operations. A full
description of Groovy semantics can be found in the [Groovy Documentation](https://groovy-lang.org/semantics.html).

## Glossary

|Term|Description|
|----|-----------|
|Domain Specific Language|A domain-specific language (DSL) is a computer language specialized to a particular application domain. Nextflow is a language for computational workflows |


FIXME

{% include links.md %}
