---
title: "Nextflow scripting"
teaching: 30
exercises: 5
questions:
- "What language are Nextflow scripts written in ? "
- "How do I assign variables ?"
- "How do I write comments ?"
- "How can I store many values together ?"
- "What are closures?"
- "How are strings evaluated in Nextflow?"
- "How can I create simple re-useabel code blocks?"
objectives:
- "Understand what language Nextflow scripts are written in."
- "Define variables in a script."
- "Create and index lists of simple values."
- "Comment Nextflow scripts."
- "Explain what a list is."
- "Explain what a string interpolation is."
- "Understand what a closure is."

keypoints:
- "Nextflow is a Domain Specific Language (DSL) implemented on top of the Groovy programming language."
- "To define a variable, simply assign a value to it e.g a = 1 ."
- "Comments use the same syntax as in the C-family programming languages: `//` or multiline `/* */`. "
- "[value1, value2, value3, ...] creates a list."
- "Lists are indexed and sliced with square brackets (e.g., list[0] and list[2..9])"
- "String interpolation (or variable interpolation, variable substitution, or variable expansion) is the process of evaluating a string literal containing one or more placeholders, yielding a result in which the placeholders are replaced with their corresponding values."
- "A closure is is a block of code that can be passed as an argument to a function."

---

Nextflow is a Domain Specific Language (DSL) implemented on top of the Groovy programming lang, which in turns is a super-set of the Java programming language. This means that Nextflow can run any Groovy and Java code. However, it is not necessary to learn Groovy to use Nextflow DSL but it can be useful in edge cases where you need more functionality than the DSL provides.

## Nextflow console

Nextflow has a console graphical interface. The  console is a REPL (read-eval-print loop) environment that allows one to quickly test part of a script or pieces of Nextflow code in an interactive manner.

It is a handy tool that allows one to evaluate fragments of Nextflow/Groovy code or fast prototype a complete pipeline script. More information can be found [here](https://www.nextflow.io/blog/2015/introducing-nextflow-console.html)

We can use the command `nextflow console` to launch the interactive console to test out out Groovy code.

~~~
nextflow console
~~~
{: .language-bash }

> ## Console Global scope
> It is worth noting that the global script context is maintained across script executions. This means that variables declared in the global script scope are not lost when the script run is complete, and they can be accessed in further executions of the same or another piece of code.
{: .callout }

## Language Basics

### Printing values

To print something is as easy as using one of the `print` or `println` methods.

~~~
println("Hello, World!")
~~~
{: .language-groovy }


The only difference between the two is that the `println` method implicitly appends a new line character to the printed string.

**Parenthesis** for function invocations are optional. Therefore also the following is a valid syntax.

~~~
println "Hello, World!"
~~~
{: .language-groovy }


## Comments

Comments use the same syntax as in the C-family programming languages. This can be confusing for people familiar with the `#` syntax for commenting in other languages.

~~~
// This is a single line comment. Everything after the // is ignored.

/*
   Comments can also
   span multiple
   lines.
 */
~~~
{: .language-groovy }


### Variables

Variables are assigned using `=` and can have any value. Variables used inside a double quoted string `""` are prefixed with a `$` to denote the variable should be interpolated, otherwise, a variable is just referenced by it’s name.

~~~
//int − This is used to represent whole numbers.
x = 1
println x

x = new java.util.Date()
println x

//float − This is used to represent 32-bit floating point numbers.
x = -3.1499392
println x

//Boolean − This represents a Boolean value which can either be true or false.
x = false
println x

//String - These are text literals which are represented in the form of chain of characters
x = "chr1"
println x
~~~
{: .language-groovy }

Local variables are defined using the def keyword:

~~~
def x = 'foo'
~~~
{: .language-groovy }

It should be always used when defining variables local to a function or a closure.

> ## Dynamic typing
> Groovy is dynamically-typed and determines its variables' data types based on their values.
Dynamically-typed languages are more flexible and can save you time and space when writing scripts. However, this can lead to issues at runtime. For example:
> ~~~
> number = 5
> numbr = (number + 15) / 2  
> ~~~
> {: .language-groovy }
> The code above should create the variable `number` with a value of 5, then change its value to 10 by adding 15 to it and dividing it by 2. However, number is misspelled at the beginning of the second line. Because Groovy does not require you to declare your variables, it creates a new variable called `numbr` and assigns it the value `number` should have. This code will compile just fine, but may produce an error later on when the script tries to do something with number assuming its value is 10.
{: .callout }


> ## Check Your Understanding
>
> What values do the variables `mass` and `age` have after each of the following statements?
> Test your answer by executing the lines.
>
> ~~~
> mass = 47.5
> age = 122
> mass = mass * 2.0
> ag = age - 20
> ~~~
> {: .language-groovy}
>
> > ## Solution
> > ~~~
> > `mass` holds a value of 47.5, `age` does not exist
> > `mass` still holds a value of 47.5, `age` holds a value of 122
> > `mass` now has a value of 95.0, `age`'s value is still 122
> > `mass` still has a value of 95.0, `age` value is still 122 as there is a typo.
> > ~~~
> > {: .output}
> {: .solution}
{: .challenge}

### Lists
A List  (also known as array) object can be defined by placing the list items in square brackets and separating items by commas `,`:

~~~
list = [10,20,30,40]
~~~
{: .language-groovy }


You can access a given item in the list with square-bracket notation or using the `get` method. These positions are numbered starting at 0, so the first element has an index of 0.

~~~
list = [10,20,30,40]
println("first element: ${list[0]}")
println("last element: ${list.get(3)}")

//Lists can also be indexed with negative indexes and ranges.
println("last element: ${list[-1]}")
println "The first three elements in the list are: ${list[0..2]}"
~~~
{: .language-groovy }

~~~
first element: 10
last element: 40
last element: 40
The first three elements in the list are: [10, 20, 30]
~~~
{: .output}

> ## Assert
> Assertions are important for checking some conditions are met.
 To perform assert add `assert` in front of variable and add `==`
 If the condition results to true, the code will continue to execute and nothing is displayed in the console and no exceptions are thrown.
{: .callout }

### List Methods

List objects implements all methods provided by the Java `java.util.List` interface plus the extension methods provided by [Groovy API](http://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/List.html).

For example, In order to get the length of the list use the `size` method:

~~~
list = [0,1,2]
println "list size is:  ${list.size()}"
~~~
{: .language-groovy }

~~~
list size is:  3
~~~
{: .output }

Listed below are a few more common list methods.

~~~
list = [1,2,3]
println list
println list + [1]
println list - [1]
println list * 2
println list.reverse()
println list.collect{ it+3 }
println list.unique().size()
println list.count(1)
println list.min()
println list.max()
println list.sum()
println list.sort()
println list.find{it%2 == 0}
println list.findAll{it%2 == 0}
~~~
{: .language-groovy }

~~~
[1, 2, 3]
[1, 2, 3, 1]
[2, 3]
[1, 2, 3, 1, 2, 3]
[3, 2, 1]
[4, 5, 6]
3
1
1
3
6
[1, 2, 3]
2
[2]
~~~
{: .output }

> ## Create List and retrieve value
> Create a list object `list` with the values 1 to 10.
> Access the fifth element in the list using with square-bracket notation or using the `get` method and
> print the results
> > ## Solution
> > ~~~
> > list = [1,2,3,4,5,6,7,8,9,10]
> > //or
> > list = 1..10
> > println "${list[4]}"
> > //or
> > println "${list.get(4)}""
> > ~~~
> > {: .language-groovy }
> > The fifth element is `5`. Remember that the array index starts at 0.
> > {: .output}
> {: .solution}
{: .challenge}


### Maps

Maps (also known as associative arrays)  are like lists that have an arbitrary type of key instead of integer. Therefore, the syntax is very much aligned. To specify the key use a colon before the value `[:]`.

~~~
emptyMap = [:]                      // an empty map
roi = [ chromosome : "chr17", start: 7640755, end: 7718054, genes: ['ATP1B2','TP53','WRAP53']]
~~~
{: .language-groovy }


Maps can be accessed in a conventional square-bracket syntax or as if the key was a property of the map
using the dor notation.

~~~
println roi['chromosome']         //Use of the square brackets.   
println roi.start                 //Use a dot notation          
println roi.get('genes')          //Use of get method
~~~
{: .language-groovy }

To add data or to modify a map, the syntax is similar to adding values to list:

~~~
roi['chromosome'] = '17'    //Use of the square brackets       
roi.chromosome = 'chr17'    //Use a dot notation          
roi.put('genome', 'hg38')   //Use of put method     
~~~
{: .language-groovy }

Map objects implements all methods provided by the `Java java.util.Map` interface plus the extension methods provided by [Groovy API](http://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Map.html).


## String interpolation


String literals can be defined enclosing them either with single-quoted `''` or double-quotes characters `""`.

Double-quoted strings can contain the value of an arbitrary variable by prefixing its name with the `$` character, or the value of any expression by using the `${expression}` syntax, similar to Bash/shell scripts:

~~~

println "The regions of interest is ${roi.chr} and contains the genes ${roi.genes.join(',')}"
println '${roi.chr}'
~~~
{: .language-groovy }


This code prints:

~~~
The regions of interest is null and contains the genes ATP1B2,TP53,WRAP53
${roi.chr}
~~~
{: .output}

Note the different use of `$` and `${..}` syntax to interpolate value expressions in a string literal.

### Slashy strings

Finally string literals can also be defined using the `/` character as delimiter. They are known as slashy strings and are useful for defining regular expressions and patterns, as there is no need to escape backslashes. As with double quote strings they allow to interpolate variables prefixed with a `$` character.

Try the following to see the difference:

~~~
x = /ATP1B2\TP53\WRAP53/
println x
~~~
{: .language-groovy }

~~~
ATP1B2\TP53\WRAP53
~~~
{: .output }

~~~
y = 'ATP1B2\TP53\WRAP53'
println y
~~~
{: .language-groovy }

Produces an error as the `\` is a special characters that we need to escape.

~~~
y = 'ATP1B2\\TP53\\WRAP53'
println y
~~~
{: .language-groovy }

~~~
ATP1B2\TP53\WRAP53
~~~
{: .output }


## Multi-line strings

A block of text that span multiple lines can be defined by delimiting it with triple single `'''` or double quotes `"""`:

~~~
text = """
    This is a multi-line
    using triple quotes.
    """
~~~
{: .language-groovy }


Finally multi-line strings can also be defined with slashy string `/`. For example:

~~~
text = /
    This is a multi-line
    slashy string.
    /
~~~
{: .language-groovy }


Like before, multi-line strings  inside double quotes `""` and slash `/` characters support variable interpolation, while single-quoted `''` multi-line strings do not.


## Closures

Closures are the swiss army knife of Nextflow/Groovy programming. In a nutshell a closure is is a block of code that can be passed as an argument to a function, it could also be defined an anonymous function.

We can assign a closure to a variable.

~~~
square = { it * it }
~~~
{: .language-groovy }


The curly brackets `{}` around the expression `it * it` tells the script interpreter to treat this expression as code. The `it` identifier is an implicit variable that represents the value that is passed to the function when it is invoked.


We can pass the function `square` as an argument to other functions or methods. Some built-in functions take a function like this as an argument. One example is the `collect` method on lists:

~~~
x = [ 1, 2, 3, 4 ]
x.collect(square)
println x
~~~
{: .language-groovy }

~~~
[ 1, 4, 9, 16 ]
~~~
{: .output }

By default, closures take a single parameter called `it`, to give it a different name use the -> syntax. For example:

~~~
square = { num -> num * num }
~~~
{: .language-groovy }

> ## Write a closure
> Write a closure to add the prefix `chr` to each element of  the list `x=[1,2,3,4,5,6]`
> > ## Solution
> > ~~~
> > prefix = {num -> "chr${num}"}
> > x = [ 1, 2, 3, 4 ].collect(prefix)
> > println x
> > ~~~
> > {: .language-groovy}
> > ~~~
> > [chr1, chr2, chr3, chr4]
> > ~~~
> > {: .output}
> {: .solution}
{: .challenge}

It’s also possible to define closures with multiple, custom-named parameters.

For example, the method `each()` when applied to a `map` can take a closure with two arguments, to which it passes the *key-value* pair for each entry in the map object. For example:

~~~
//closure
printMap = { a, b -> println "$a with value $b" }

//map object
values = [ "chromsome" : "chr17", "start" : 1, "end" : 83257441 ]
values.each(printMap)
~~~
{: .language-groovy }


~~~
chromsome with value chr17
start with value 1
end with value 83257441
~~~
{: .output }


A closure has two other important features. First, it can access and modify variables in the scope where it is defined.

Second, a closure can be defined in an anonymous manner, meaning that it is not given a name, and is defined in the place where it needs to be used.

As an example showing both these features, see the following code fragment:

~~~
result = 0  
values = ["China": 1 , "India" : 2, "USA" : 3]   
values.keySet().each({ result += values[it] })    
println result
~~~
{: .language-groovy }


* Define a global variable.
* Define a map object.
* Invoke the each method passing closure object which modifies the result variable.

Learn more about closures in the [Groovy documentation](http://groovy-lang.org/closures.html).




# Additional Material

## Conditional Execution

## If statement

One of the most important features of any programming language is the ability to execute different code under different conditions. The simplest way to do this is to use the if construct:

The if statement uses the same syntax common other programming lang such Java, C, JavaScript, etc.

~~~
if( < boolean expression > ) {
    // true branch
}
else {
    // false branch
}
~~~
{: .language-groovy }


The else branch is optional. Also curly brackets are optional when the branch define just a single statement.

~~~
x = 12
if( x > 10 )
    println "$x is greater the 10"
~~~
{: .language-groovy }


null, empty strings and empty collections are evaluated to false.
Therefore a statement like:

~~~
list = [1,2,3]
if( list != null && list.size() > 0 ) {
  println list
}
else {
  println 'The list is empty'
}
~~~
{: .language-groovy }


Can be written as:

~~~
if( list )
    println list
else
    println 'The list is empty'
~~~
{: .language-groovy }


In some cases can be useful to replace `if` statement with a ternary expression aka conditional expression. For example:

~~~
println list ? list : 'The list is empty'
~~~
{: .language-groovy }


The previous statement can be further simplified using the Elvis operator `?:` as shown below:

~~~
println list ?: 'The list is empty'
~~~
{: .language-groovy }


## For statement

The classical for loop syntax is supported as shown here:

~~~
for (int i = 0; i <3; i++) {
   println("Hello World $i")
}
~~~
{: .language-groovy }


Iteration over list objects is also possible using the syntax below:

~~~
list = ['a','b','c']

for( String elem : list ) {
  println elem
}
~~~
{: .language-groovy }


## Functions

It is possible to define a custom function into a script, as shown here:

~~~
int fib(int n) {
    return n < 2 ? 1 : fib(n-1) + fib(n-2)
}

assert fib(10)==89
~~~
{: .language-groovy }


A function can take multiple arguments separating them with a comma. The return keyword can be omitted and the function implicitly returns the value of the last evaluated expression. Also explicit types can be omitted (thought not recommended):

~~~
def fact( n ) {
  n > 1 ? n * fact(n-1) : 1
}


assert fact(5) == 120
~~~
{: .language-groovy }



## More resources

The complete Groovy language documentation is available at this [link](http://groovy-lang.org/documentation.html#languagespecification).

A great resource to master Apache Groovy syntax is Groovy in [Action](https://www.manning.com/books/groovy-in-action-second-edition).
