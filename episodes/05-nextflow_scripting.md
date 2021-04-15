---
title: "Nextflow scripting"
teaching: 0
exercises: 0
questions:
- "What language are Nextflow scripts written in?"
- "How do you assign variables?"
- "How do write comments?"
- "How do you control code execution?"
- "How do you loop through lists?"
- "How do write reuseable code blocks?"
- "What are closures?"
- "How are strings evaluated in Nextflow?"
objectives:
- "Understand a Nextflow process."
- "Comment Nextflow scripts."
- "Control the execution using conditional statements."

keypoints:
- "Nextflow is a Domain Specific Language (DSL) implemented on top of the Groovy programming language."
- "To define a variable, simply assign a value to it e.g a = 1 ."

---

Nextflow is a Domain Specific Language (DSL) implemented on top of the Groovy programming lang, which in turns is a super-set of the Java programming language. This means that Nextflow can run any Groovy and Java code. However, it is not necessary to learn Groovy to use Nextflow DSL but it can be useful.

## Language Basics

### Printing values

To print something is as easy as using one of the `print` or `println` methods.

~~~
println("Hello, World!")
~~~
{: .source}

The only difference between the two is that the `println` method implicitly appends a new line character to the printed string.

**parenthesis** for function invocations are optional. Therefore also the following is a valid syntax.

~~~
println "Hello, World!"
~~~
{: .source}

## Comments

Comments use the same syntax as in the C-family programming languages. This can be confusing for people familiar with the `#` syntax for commenting in other languages.

~~~
// comment a single config file

/*
   a comment spanning
   multiple lines
 */
~~~
{: .source}

### Variables


To define a variable, simply assign a value to it using  `=`:

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
x = "Hi"
println x
~~~
{: .source}

Local variables are defined using the `def` keyword:

~~~
def x = 'foo'
~~~
{: .source}

It should be always used when defining variables local to a function or a closure.

### Lists

A List object can be defined by placing the list items in square brackets `[]` and separating items by commas `,`:

~~~
list = [10,20,30,40]
~~~
{: .source}

You can access a given item in the list with square-bracket notation (indexes start at `0`) or using the `get` method:

~~~
assert list[0] == 10
assert list[0] == list.get(0)
~~~
{: .source}

> ## Assert
> Assertions are important for checking some conditions are met.
> To perform assert add `assert` in front of variable and add `==`
> If the condition results to true, the code will continue to execute and nothing is displayed in the console and no exceptions are thrown.
> {: .callout}

In order to get the length of the list use the `size` method:


~~~
assert list.size() == 4
Lists can also be indexed with negative indexes and reversed ranges.

list = [0,1,2]
assert list[-1] == 2
assert list[-1..0] == list.reverse()
~~~
{: .source}

List objects implements all methods provided by the Java java.util.List interface plus the extension methods provided by [Groovy API](http://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/List.html).


~~~
assert [1,2,3] << 1 == [1,2,3,1]
assert [1,2,3] + [1] == [1,2,3,1]
assert [1,2,3,1] - [1] == [2,3]
assert [1,2,3] * 2 == [1,2,3,1,2,3]
assert [1,[2,3]].flatten() == [1,2,3]
assert [1,2,3].reverse() == [3,2,1]
assert [1,2,3].collect{ it+3 } == [4,5,6]
assert [1,2,3,1].unique().size() == 3
assert [1,2,3,1].count(1) == 2
assert [1,2,3,4].min() == 1
assert [1,2,3,4].max() == 4
assert [1,2,3,4].sum() == 10
assert [4,2,1,3].sort() == [1,2,3,4]
assert [4,2,1,3].find{it%2 == 0} == 4
assert [4,2,1,3].findAll{it%2 == 0} == [4,2]
~~~
{: .source}



### Maps


Maps are like lists that have an arbitrary type of key instead of integer. Therefore, the syntax is very much aligned. To specify the key use a colon before the value `[:]`.

~~~
map = [a:0, b:1, c:2]
~~~
{: .source}

Maps can be accessed in a conventional square-bracket syntax or as if the key was a property of the map.

~~~
assert map['a'] == 0        
assert map.b == 1           
assert map.get('c') == 2
~~~

* Use of the square brackets.
* Use a dot notation.
* Use of get method.

To add data or to modify a map, the syntax is similar to adding values to list:

~~~
map['a'] = 'x'           
map.b = 'y'              
map.put('c', 'z')        
assert map == [a:'x', b:'y', c:'z']
~~~
{: .source}

* Use of the square brackets.
* Use a dot notation.
* Use of get method.


Map objects implements all methods provided by the Java java.util.Map interface plus the extension methods provided by [Groovy API](http://docs.groovy-lang.org/latest/html/groovy-jdk/java/util/Map.html).


## String interpolation

String literals can be defined enclosing them either with single-quoted `''` or double-quotes characters `""`.

Double-quoted strings can contain the value of an arbitrary variable by prefixing its name with the `$` character, or the value of any expression by using the `${expression}` syntax, similar to Bash/shell scripts:

~~~
foxtype = 'quick'
foxcolor = ['b', 'r', 'o', 'w', 'n']
println "The $foxtype ${foxcolor.join()} fox"

x = 'Hello'
println '$x + $y'
~~~
{: .source}

This code prints:

~~~
The quick brown fox
$x + $y
~~~
{: .ouput}

Note the different use of `$` and `${..}` syntax to interpolate value expressions in a string literal.
Finally string literals can also be defined using the `/` character as delimiter. They are known as slashy strings and are useful for defining regular expressions and patterns, as there is no need to escape backslashes. As with double quote strings they allow to interpolate variables prefixed with a `$` character.

Try the following to see the difference:

~~~
x = /tic\tac\toe/
y = 'tic\tac\toe'

println x
println y
~~~
{: .source}

it prints:

~~~
tic\tac\toe
tic    ac    oe
~~~
{: .source}

## Multi-line strings

A block of text that span multiple lines can be defined by delimiting it with triple single `'''` or double quotes `"""`:

~~~
text = """
    Hello there James
    how are you today?
    """
~~~
{: .source}

Finally multi-line strings can also be defined with slashy string `/`. For example:

~~~
text = /
    This is a multi-line
    slashy string!
    It's cool, isn't it?!
    /
~~~
{: .source}

Like before, multi-line strings inside double quotes and slash characters support variable interpolation, while single-quoted multi-line strings do not.

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
{: .source}

The else branch is optional. Also curly brackets are optional when the branch define just a single statement.

~~~
x = 1
if( x > 10 )
    println 'Hello'
~~~
{: .source}

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
{: .source}

Can be written as:

~~~
if( list )
    println list
else
    println 'The list is empty'
~~~
{: .source}

See the Groovy-Truth for details.

In some cases can be useful to replace if statement with a ternary expression aka conditional expression. For example:

~~~
println list ? list : 'The list is empty'
~~~
{: .source}

The previous statement can be further simplified using the Elvis operator `?:` as shown below:

~~~
println list ?: 'The list is empty'
~~~
{: .source}

## For statement

The classical for loop syntax is supported as shown here:

~~~
for (int i = 0; i <3; i++) {
   println("Hello World $i")
}
~~~
{: .source}

Iteration over list objects is also possible using the syntax below:

~~~
list = ['a','b','c']

for( String elem : list ) {
  println elem
}
~~~
{: .source}

## Functions

It is possible to define a custom function into a script, as shown here:

~~~
int fib(int n) {
    return n < 2 ? 1 : fib(n-1) + fib(n-2)
}


assert fib(10)==89
~~~
{: .source}

A function can take multiple arguments separating them with a comma. The return keyword can be omitted and the function implicitly returns the value of the last evaluated expression. Also explicit types can be omitted (thought not recommended):

~~~
def fact( n ) {
  n > 1 ? n * fact(n-1) : 1
}


assert fact(5) == 120
~~~
{: .source}

## Closures

Closures are the swiss army knife of Nextflow/Groovy programming. In a nutshell a closure is is a block of code that can be passed as an argument to a function, it could also be defined an anonymous function.

More formally, a closure allows the definition of functions as first class objects.

~~~
square = { it * it }
~~~
{: .source}

The curly brackets `{}` around the expression `it * it` tells the script interpreter to treat this expression as code. The `it` identifier is an implicit variable that represents the value that is passed to the function when it is invoked.

Once compiled the function object is assigned to the variable square as any other variable assignments shown previously. To invoke the closure execution use the special method `call` or just use the round parentheses to specify the closure parameter(s). For example:

~~~
assert square.call(5) == 25
assert square(9) == 81
~~~
{: .source}

This is not very interesting until we find that we can pass the function square as an argument to other functions or methods. Some built-in functions take a function like this as an argument. One example is the collect method on lists:

~~~
x = [ 1, 2, 3, 4 ].collect(square)
println x
~~~
{: .source}

It prints:

~~~
[ 1, 4, 9, 16 ]
~~~
{: .source}

By default, closures take a single parameter called `it`, to give it a different name use the -> syntax. For example:

~~~
square = { num -> num * num }
~~~

It’s also possible to define closures with multiple, custom-named parameters.

For example, the method `each()` when applied to a `map` can take a closure with two arguments, to which it passes the *key-value* pair for each entry in the map object. For example:

~~~
printMap = { a, b -> println "$a with value $b" }
values = [ "Yue" : "Wu", "Mark" : "Williams", "Sudha" : "Kumari" ]
values.each(printMap)
~~~
{: .source}

It prints:

~~~
Yue with value Wu
Mark with value Williams
Sudha with value Kumari
~~~
{: .source}

A closure has two other important features. First, it can access and modify variables in the scope where it is defined.

Second, a closure can be defined in an anonymous manner, meaning that it is not given a name, and is defined in the place where it needs to be used.

As an example showing both these features, see the following code fragment:

~~~
result = 0  
values = ["China": 1 , "India" : 2, "USA" : 3]   
values.keySet().each({ result += values[it] })    
println result
~~~
{: .source}

* Define a global variable.
* Define a map object.
* Invoke the each method passing closure object which modifies the result variable.

Learn more about closures in the [Groovy documentation](http://groovy-lang.org/closures.html).

## More resources

The complete Groovy language documentation is available at this [link](http://groovy-lang.org/documentation.html#languagespecification).

A great resource to master Apache Groovy syntax is Groovy in [Action](https://www.manning.com/books/groovy-in-action-second-edition).
