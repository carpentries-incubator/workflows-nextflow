---
title: Nextflow scripting
teaching: 30
exercises: 5
---

::::::::::::::::::::::::::::::::::::::: objectives

- Understand what language Nextflow scripts are written in.
- Define variables in a script.
- Create lists of simple values.
- Comment Nextflow scripts.
- Explain what a list is.
- Explain what string interpolation is.
- Understand what a closure is.

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::: questions

- What language are Nextflow scripts written in?
- How do I store values in a Nextflow script?
- How do I write comments in a Nextflow script?
- How can I store and retrieve multiple values?
- How are strings evaluated in Nextflow?
- How can I create simple re-useable code blocks?

::::::::::::::::::::::::::::::::::::::::::::::::::

Nextflow is a Domain Specific Language (DSL) implemented on top of the Groovy programming language, which in turn is a super-set of the Java programming language. This means that Nextflow can run any Groovy and Java code. It is not necessary to learn Groovy to use Nextflow DSL but it can be useful in edge cases where you need more functionality than the DSL provides.

::::::::::::::::::::::::::::::::::::::::  callout

## Nextflow console

Nextflow has a console graphical interface. The  console is a REPL (read-eval-print loop) environment that allows a user to quickly test part of a script or pieces of Nextflow code in an interactive manner.

It is a handy tool that allows a user to evaluate fragments of Nextflow/Groovy code or fast prototype a complete pipeline script. More information can be found [here](https://www.nextflow.io/blog/2015/introducing-nextflow-console.html)

We can use the command `nextflow console` to launch the interactive console to test out out Groovy code.

```bash 
nextflow console
```

### Console global scope

It is worth noting that the global script context is maintained across script executions. This means that variables declared in the global script scope are not lost when the script run is complete, and they can be accessed in further executions of the same or another piece of code.


::::::::::::::::::::::::::::::::::::::::::::::::::

## Language Basics

### Printing values

To print something is as easy as using the `println` method (`println` is a compression of "print line") and passing the text to print in quotes.
The text is referred to as a `string` as in a string of characters.

```groovy 
println("Hello, World!")
```

```output 
Hello, World!
```

**Parenthesis** for function invocations are optional. Therefore also the following is a valid syntax.

```groovy 
println "Hello, World!"
```

```output 
Hello, World!
```

## Methods

`println` is a example of a Groovy method. A method is just a block of code which only runs when it is called.
You can pass data, known as parameters, into a method using the method name followed by brackets `()`.
Methods are used to perform certain actions, and they are also known as functions.
Methods enable us to reuse code: define the code once, and use it many times.

## Comments

When we write any code it is useful to document it using comments.
In Nextflow comments use the same syntax as in the C-family programming languages.
This can be confusing for people familiar with the `#` syntax for commenting in other languages.

```groovy 
// This is a single line comment. Everything after the // is ignored.

/*
   Comments can also
   span multiple
   lines.
 */
```

## Variables

In any programming language, you need to use variables to store different types of information. A variable is a pointer to a space in the computer's memory that stores the value associated with it.

Variables are assigned using `=` and can have any value. Groovy is dynamically-typed which means the variable's data type is based on its value. For example, setting `x = 1` means `x` is an integer number, but if it is later set to `x = "hello"` then it becomes a String.

::::::::::::::::::::::::::::::::::::::::  callout

## Variable scope

When we create a variable using the `x = 1` syntax we can access, (`scope`), it anywhere (`globally`) in the script. A variable declared in this fashion is sometimes called a public variable.

We can also define variables with a data `type` e.g. `String x="Hello"` or with the `def` keyword `def x=1`. This effects the accessibility (`scope`) of the variable.
This is called lexical scoping (sometimes known as static scoping) that sets the scope  of a variable so that it may only be accessed from within the block of code in which it is defined. A variable declared in this fashion is sometimes called a private variable.


::::::::::::::::::::::::::::::::::::::::::::::::::

### Types of Data

Groovy knows various types of data. four common ones are:

- `String` − These are text literals which are represented in the form of chain of characters enclosed in quotes. For example `"Hello World"`.
- `int` − This is used to represent whole numbers. An example is `1234`.
- `Boolean` − This represents a Boolean value which can either be `true` or `false`.
- `float` - This is used to represent floating point number `12.34` .

A more complete list can be found [here](https://www.tutorialspoint.com/groovy/groovy_data_types.htm)

In the example below, variable `my_var` has an integer value of `1`:

```groovy 
//int − This is used to represent whole numbers.
my_var = 1
```

To create a variable with a floating point value, we can execute:

```groovy 
//float − This is used to represent floating point numbers.
my_var = 3.1499392
```

To create a Boolean value we assign the value `true` or `false`.  
\**Note:* Do not enclose a Boolean value in quotes or they will be interpreted as a string.

```groovy 
//Boolean − This represents a Boolean value which can either be true or false.
my_var = false
```

And to create a string, we add single or double quotes around some text.

For example:

```groovy 
//String - These are text literals which are represented in the form of chain of characters
my_var = "chr1"
```

### Multi-line strings

A block of text that span multiple lines can be defined by delimiting it with triple single `'''` or double quotes `"""`:

```groovy 
text = """
    This is a multi-line string
    using triple quotes.
    """
```

To display the value of a variable to the screen in Groovy, we can use the `println` method passing the variable name are a parameter.

```groovy 
x = 1
println(x)
```

```output 
1
```

### Slashy strings

Strings can also be defined using the forward slash `/` character as delimiter. They are known as `slashy strings` and are useful for defining regular expressions and patterns, as there is no need to escape backslashes e.g `\n` specifies a new line. As with double quote strings they allow to interpolate variables prefixed with a `$` character.

Try the following to see the difference:

```groovy 
x = /ATP1B2\TP53\WRAP53/
println(x)
```

```output 
ATP1B2\TP53\WRAP53
```

```groovy 
y = 'ATP1B2\TP53\WRAP53'
println(y)
```

Produces an error as the `\` is a special characters that we need to escape.

```groovy 
// use \ to escape
y = 'ATP1B2\\TP53\\WRAP53'
println(y)
```

```output 
ATP1B2\TP53\WRAP53
```

### String interpolation

To use a variable inside a single or multi-line double quoted string `""`  prefix the variable name with a `$` to show it should be interpolated:

```groovy 
chr = "1"
println("processing chromosome $chr")
```

```output 
processing chromosome 1
```

**Note:** Variable names inside single quoted strings do not support String interpolation.

```groovy 
chr = "1"
println('processing chromosome $chr')
```

```output 
processing chromosome $chr
```

## Lists

To store multiple values in a variable we can use a List.
A List  (also known as array) object can be defined by placing the list items in square brackets and separating items by commas `,`:

```groovy 
kmers = [11,21,27,31]
```

You can access a given item in the list with square-bracket notation `[]`. These positions are numbered starting at 0, so the first element has an index of 0.

```groovy 
kmers = [11,21,27,31]
println(kmers[0])

```

```output
11
```

We can use negative numbers as indices in Groovy. They count from the end of the list rather than the front: the index `-1` gives us the last element in the list, `-2` the second to last, and so on. Because of this, `kmers[3]` and `kmers[-1]` point to the same element in our example list.

```groovy 
kmers = [11,21,27,31]
//Lists can also be indexed with negative indexes
println(kmers[3])
println(kmers[-1])
```

```output
31
31
```

Lists can also be indexed using a range. A range is a quick way of declaring a list of consecutive sequential numbers.
To define a range use `<num1>..<num2>` notation.

```groovy 
kmers = [11,21,27,31]
// The first three elements using a range.
println(kmer[0..2])
```

```output
[11, 21, 27]
```

### String interpolation of list elements

To use an expression like `kmer[0..2]` inside a double quoted String `""` we use the `${expression}` syntax, similar to Bash shell scripts.

For example, the expression below without the `{}`""

```groovy 
kmers = [11,21,27,31]
println("The first three elements in the Lists are. $kmers[0..2]")
```

would output.

```output
The first three elements in the Lists are. [11, 21, 27, 31][0..2]
```

We need to enclose the `kmers[0..2]` expression inside `${}` as below to get the correct output.

```groovy 
kmers = [11,21,27,31]
println("The first three elements in the Lists are. ${kmers[0..2]}")
```

```output
The first three elements in the Lists are. [11, 21, 27]
```

### List methods

Lists have a number of useful methods that can perform operations on their contents. See more [here](https://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/List.html). When using a method on a type of object you need prefix the method with the variable name.

For example, in order to get the length of the list use the list `size` method:

```groovy 
mylist = [0,1,2]

println(mylist.size())

//inside a string need we need to use the ${} syntax
println("list size is:  ${mylist.size()}")
```

```output 
3
list size is:  3
```

We can use the `get` method items to retrieve items in a list.

```groovy 
mylist = [0,1,2]
println mylist.get(1)
```

```output 
1
```

Listed below are a few more common list methods and their output on a simple example.

```groovy 
mylist = [1,2,3]
println mylist
println mylist + [1]
println mylist - [1]
println mylist * 2
println mylist.reverse()
println mylist.collect{ it+3 }
println mylist.unique().size()
println mylist.count(1)
println mylist.min()
println mylist.max()
println mylist.sum()
println mylist.sort()
println mylist.find{it%2 == 0}
println mylist.findAll{it%2 == 0}
```

```output 
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
```

:::::::::::::::::::::::::::::::::::::::  challenge

## Create List and retrieve value

Create a list object `list` with the values 1 to 10.
Access the fifth element in the list using with square-bracket notation or using the `get` method and
print the results

:::::::::::::::  solution

## Solution

```groovy 
list = [1,2,3,4,5,6,7,8,9,10]
//or
list = 1..10
println("${list[4]}")
//or
println("${list.get(4)}")
```

The fifth element is `5`. Remember that the array index starts at 0.



:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

## Maps

It can difficult to remember the index of a value in a list, so we can use Groovy Maps (also known as associative arrays) that have an arbitrary type of key instead of an integer value. The syntax is very similar to Lists. To specify the key use a colon before the value `[key:value]`. Multiple values are separated by a comma. *Note:* the key value is not enclosed in quotes.

```groovy 
roi = [ chromosome : "chr17", start: 7640755, end: 7718054, genes: ['ATP1B2','TP53','WRAP53']]
```

Maps can be accessed in a conventional square-bracket syntax or as if the key was a property of the map or using the dot notation. *Note:* When retrieving a value the key value is  enclosed in quotes.

```groovy 
//Use of the square brackets.
println(roi['chromosome'])

//Use a dot notation            
println(roi.start)

//Use of get method                      
println(roi.get('genes'))          
```

To add data or to modify a map, the syntax is similar to adding values to list:

```groovy 
//Use of the square brackets
roi['chromosome'] = '17'

//Use a dot notation        
roi.chromosome = 'chr17'  

//Use of put method              
roi.put('genome', 'hg38')  
```

More information about maps can be found in the [Groovy API](https://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Map.html).

## Closures

Closures are the swiss army knife of Nextflow/Groovy programming. In a nutshell a closure is a block of code that can be passed as an argument to a function. This can be useful to create a re-usable function.

We can assign a closure to a variable in same way as a value using the `=`.

```groovy 
square = { it * it }
```

The curly brackets `{}` around the expression `it * it` tells the script interpreter to treat this expression as code. `it` is an implicit variable that is provided in closures. It's available when the closure doesn't have an explicitly declared parameter and represents the value that is passed to the function when it is invoked.

We can pass the function `square` as an argument to other functions or methods. Some built-in functions take a function like this as an argument. One example is the `collect` method on lists that iterates through each element of the list transforming it into a new value using the closure:

```groovy 
square = { it * it }
x = [ 1, 2, 3, 4 ]
y = x.collect(square)
println y
```

```output 
[ 1, 4, 9, 16 ]
```

A closure can also be defined in an anonymous manner, meaning that it is not given a name, and is defined in the place where it needs to be used.

```groovy 
x = [ 1, 2, 3, 4 ]
y = x.collect({ it * it })
println("x is $x")
println("y is $y")
```

```output 
x is [1, 2, 3, 4]
y is [1, 4, 9, 16]
```

### Closure parameters

By default, closures take a single parameter called `it`. To define a different name use the
` variable ->` syntax.

For example:

```groovy 
square = { num -> num * num }
```

In the above example the variable `num` is assigned as the closure input parameter instead of `it`.

:::::::::::::::::::::::::::::::::::::::  challenge

## Write a closure

Write a closure to add the prefix `chr` to each element of  the list `x=[1,2,3,4,5,6]`

:::::::::::::::  solution

## Solution

```groovy
prefix = { "chr${it}"}
x = [ 1,2,3,4,5,6 ].collect(prefix)
println x
```

```output
[chr1, chr2, chr3, chr4, chr5, chr6]
```

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

### Multiple map parameters

It's also possible to define closures with multiple, custom-named parameters using the `->` syntax. This separate the custom-named parameters by a comma before the `->` operator.

For example:

```groovy 
tp53 = [chromosome: "chr17",start:7661779 ,end:7687538, genome:'GRCh38', gene: "TP53"]
//perform subtraction of end and start coordinates
region_length = {start,end -> end-start }
tp53.length = region_length(tp53.start,tp53.end)
println(tp53)
```

Would add the region `length` to the map `tp53`, calculated as `end - start`.

```output 
[chromosome:chr17, start:7661779, end:7687538, genome:GRCh38, gene:TP53, length:25759]
```

For another example, the method `each()` when applied to a `map` can take a closure with two arguments, to which it passes the *key-value* pair for each entry in the map object:

```groovy 
//closure with two parameters
printMap = { a, b -> println "$a with value $b" }

//map object
my_map = [ chromosome : "chr17", start : 1, end : 83257441 ]

//each iterates through each element
my_map.each(printMap)
```

```output 
chromosome with value chr17
start with value 1
end with value 83257441
```

Learn more about closures in the [Groovy documentation](https://groovy-lang.org/closures.html).

## Additional Material

## Conditional Execution

## If statement

One of the most important features of any programming language is the ability to execute different code under different conditions. The simplest way to do this is to use the if construct.

The if statement uses the syntax common to other programming languages such Java, C, JavaScript, etc.

```groovy 
if( < boolean expression > ) {
    // true branch
}
else {
    // false branch
}
```

The else branch is optional. Curly brackets are optional when the branch defines just a single statement.

```groovy 
x = 12
if( x > 10 )
    println "$x is greater than 10"
```

*null*, *empty strings* and *empty collections* are evaluated to false.
Therefore a statement like:

```groovy 
list = [1,2,3]
if( list != null && list.size() > 0 ) {
  println list
}
else {
  println 'The list is empty'
}
```

Can be written as:

```groovy 
if( list )
    println list
else
    println 'The list is empty'
```

In some cases can be useful to replace `if` statement with a ternary expression, also known as a conditional expression. For example:

```groovy 
println list ? list : 'The list is empty'
```

The previous statement can be further simplified using the Elvis operator `?:` as shown below:

```groovy 
println list ?: 'The list is empty'
```

## For statement

The classical for loop syntax is supported as shown here:

```groovy 
for (int i = 0; i <3; i++) {
   println("Hello World $i")
}
```

Iteration over list objects is also possible using the syntax below:

```groovy 
list = ['a','b','c']

for( String elem : list ) {
  println elem
}
```

## Functions

It is possible to define a custom function into a script, as shown here:

```groovy 
int fib(int n) {
    return n < 2 ? 1 : fib(n-1) + fib(n-2)
}

println (fib(10)) // prints 89
```

- A function can take multiple arguments separated by commas.
- The `return` keyword can be omitted and the function implicitly returns the value of the last evaluated expression. (Not recommended)
- Explicit types can be omitted. (Not recommended):

```groovy 
def fact( n ) {
  n > 1 ? n * fact(n-1) : 1
}

println (fact(5)) // prints 120
```

## More resources

The complete Groovy language documentation is available at this [link](https://groovy-lang.org/documentation.html#languagespecification).

A great resource to master Apache Groovy syntax is Groovy in [Action](https://www.manning.com/books/groovy-in-action-second-edition).

:::::::::::::::::::::::::::::::::::::::: keypoints

- Nextflow is a Domain Specific Language (DSL) implemented on top of the Groovy programming language.
- To define a variable, assign a value to it e.g., `a = 1` .
- Comments use the same syntax as in the C-family programming languages: `//` or multiline `/* */`.
- Multiple values can be stored in lists \[value1, value2, value3, ...\] or maps \[chromosome: 1, start :1\].
- Lists are indexed and sliced with square brackets (e.g., list\[0\] and list\[2..9\])
- String interpolation (variable interpolation, variable substitution, or variable expansion) is the process of evaluating a string literal containing one or more placeholders, yielding a result in which the placeholders are replaced with their corresponding values.
- A closure is an expression (block of code) encased in `{}` e.g. `{ it * it }`.

::::::::::::::::::::::::::::::::::::::::::::::::::


