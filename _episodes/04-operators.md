---
title: "Operators"
teaching: 0
exercises: 0
questions:
- "How can connect channels?"
objectives:
- "You can connect channels using operators"
keypoints:
- "Nextflow *operators* are methods that allow you to connect channels to each other or to transform values emitted by a channel applying some user provided rules"
---

# Operators

Nextflow *operators* are methods that allow you to connect channels to each other or to transform values emitted by a channel applying some user provided rules.

* Built-in functions applied to channels
* Transform channels content
* Can be used also to filter, fork and combine channels

## Basic example

~~~
nums = Channel.from(1,2,3,4)         
square = nums.map { it -> it * it }  
square.view()  
~~~
{: .source}

![Nextflow Channel map](../fig/channel-map.png)

* Create a queue channel emitting four values.
* Create a new channels transforming each number in it’s square.
* Print the channel content.

Operators can be chained to implement custom behaviors:

~~~
Channel.from(1,2,3,4)
        .map { it -> it * it }
        .view()

~~~
{: .source}


Operators can be separated in to five groups:

* Filtering operators
* Transforming operators
* Splitting operators
* Combining operators
* Forking operators
* Maths operators

 ## Basic operators

### view

The `view` operator prints the items emitted by a channel to the console standard output appending a *new line* character to each of them. For example:
~~~
Channel
      .from('foo', 'bar', 'baz')
      .view()
~~~
{: .source}

It prints:

~~~
foo
bar
baz
~~~
{: .source}

An optional *closure* `{}` parameter can be specified to customize how items are printed. For example:

> # Closures
> Briefly, a closure is a block of code that can be passed as an argument to a function. Thus, you can define a chunk of code and then pass it around as if it were a string or an integer.
{: .callout}


~~~
Channel
      .from('foo', 'bar', 'baz')
      .view({ "- $it" })
~~~
{: .source}

It prints:

~~~
- foo
- bar
- baz
~~~
{: .source}

### map

The `map` operator applies a function of your choosing to every item emitted by a channel, and returns the items so obtained as a new channel. The function applied is called the mapping function and is expressed with a closure as shown in the example below:

~~~
Channel
    .from( 'hello', 'world' )
    .map ({ it -> it.reverse() })
    .view()
~~~
{: .source}

A map can associate to each element a generic tuple containing any data as needed.

~~~
Channel
    .from( 'hello', 'world' )
    .map ({ word -> [word, word.size()] })
    .view ({ word, len -> "$word contains $len letters" })
~~~
{: .source}

> ## fromPath
>
> Use `fromPath` to create a channel emitting the fastq files matching the pattern `data/ggal/*.fq`, then chain with a map to return a pair containing the file name and the path itself. Finally print the resulting channel.
>
> > ## Solution
> >
> > ~~~
> > Channel
> >   .fromPath( 'data/ggal/*.fq' )
> >   .map ({file -> [ file.name, file ]})
> >   .view({name, file -> "> file: $name"})
> > ~~~    
> {: .solution}
{: .challenge}

### into
The `into` operator connects a source channel to two or more target channels in such a way the values emitted by the source channel are copied to the target channels. For example:


~~~
Channel
     .from( 'a', 'b', 'c' )
     .into{ foo; bar }

foo.view{ "Foo emits: " + it }
bar.view{ "Bar emits: " + it }
~~~
{: .source}

> ## channel names separator
> > Note the use in this example of curly brackets and the ; as channel names separator. This is needed because the actual parameter of into is a closure which defines the target channels to which the source one is connected.
> >

### mix

The `mix` operator combines the items emitted by two (or more) channels into a single channel.
~~~
c1 = Channel.from( 1,2,3 )
c2 = Channel.from( 'a','b' )
c3 = Channel.from( 'z' )

c1 .mix(c2,c3).view()
~~~
{: .source}

~~~
1
2
a
3
b
z
~~~
{: .output}

The items in the resulting channel have the same order as in respective original channel, however there’s no guarantee that the element of the second channel are append after the elements of the first. Indeed in the above example the element a has been printed before 3.


###  flatten
The `flatten` operator transforms a channel in such a way that every *tuple* is flattened so that each single entry is emitted as a sole element by the resulting channel.

~~~
foo = [1,2,3]
bar = [4, 5, 6]

Channel
    .from(foo, bar)
    .flatten()
    .view()

~~~
{: .source}
The above snippet prints:
~~~
1
2
3
4
5
6
~~~
{: .output}

### collect

The `collect` operator collects all the items emitted by a channel to a list and return the resulting object as a sole emission.

~~~
Channel
    .from( 1, 2, 3, 4 )
    .collect()
    .view()
~~~
{: .source}


It prints a single value:

The result of the collect operator is a value channel.

~~~
[1,2,3,4]
~~~
{: .output}

### groupTuple

The `groupTuple` operator collects *tuples* (or lists) of values emitted by the source channel grouping together the elements that share the same key. Finally it emits a new tuple object for each distinct key collected.

Try the following example:

~~~
Channel
     .from( [1,'A'], [1,'B'], [2,'C'], [3, 'B'], [1,'C'], [2, 'A'], [3, 'D'] )
     .groupTuple()
     .view()
~~~~     
{: .source}

It shows:
~~~
[1, [A, B, C]]
[2, [C, A]]
[3, [B, D]]
~~~

This operator is useful to process altogether all elements for which there’s a common property or a grouping key.

Exercise


> ## Challenge Title
>
> Use `fromPath` to create a channel emitting the fastq files matching the pattern `data/ggal/*.fq`, then use a map to associate to each file the name prefix. Finally group together all files having the same common prefix.
>
> > ## Solution
> >
> > ~~~
> > Channel.fromPath('data/ggal/*.fq')
> >     .map { file -> [ file.name.split('_')[0], file ] }
> >     .groupTuple()
> >     .view()
> > ~~~
> {: .solution}
{: .challenge}

### join

The `join` operator creates a channel that joins together the items emitted by two channels for which exits a matching key. The key is defined, by default, as the first element in each item emitted.

~~~
left = Channel.from(['X', 1], ['Y', 2], ['Z', 3], ['P', 7])
right= Channel.from(['Z', 6], ['Y', 5], ['X', 4])
left.join(right).view()
~~~
{: .source}

The resulting channel emits:
~~~
[Z, 3, 6]
[Y, 2, 5]
[X, 1, 4]
~~~

### branch

The `branch` operator allows you to forward the items emitted by a source channel to one or more output channels, choosing one out of them at a time.

The selection criteria is defined by specifying a closure that provides one or more boolean expression, each of which is identified by a unique label. On the first expression that evaluates to a true value, the current item is bound to a named channel as the label identifier. For example:

~~~
Channel
    .from(1,2,3,40,50)
    .branch ({
        small: it < 10
        large: it > 10
    })
    .set({ result })

 result.small.view({ "$it is small" })
 result.large.view({ "$it is large" })
~~~
{: .source}

> ## multi-channel object
> The branch operator returns a multi-channel object i.e. a variable that holds more than one channel object.

## More resources

Check the operators [documentation](https://www.nextflow.io/docs/latest/operator.html) on Nextflow web site.

{: .output}
{% include links.md %}
