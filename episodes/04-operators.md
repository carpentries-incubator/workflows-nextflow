---
title: "Operators"
teaching: 30
exercises: 10
questions:
- "What are Nextflow operators?"
- "How do you perform operations such as filtering on channels?"
- "What are the different kinds of operators?"
- "How do you combine operators?"
- "How do I process a csv file using operators?"
objectives:
- "Describe what Nextflow operators are."
- "Modify the contents/elements of a channel"
- "Perform filtering and combining operations on a channel. "
- "Use the `splitCsv` operator to process a csv file."

keypoints:
- "Nextflow *operators* are methods that allow you to modify Nextflow channels either by splitting or combining channels, or to transform elements within a channel applying function"
- "An operator is method that transforms a channel into a new one by applying a function to each element."
- "You can connect channels using operators"
- "Operators can be separated in to  groups:; filtering , transforming , splitting , combining , forking and Maths operators"
- "You can split csv file using the splitCsv operator."
---

# Operators

 *Operators* is how Nextflow allow you to modify the contents of a channel.

 Operators can be separated in to several groups:


 * Filtering operators
 * Transforming operators
 * Splitting operators
 * Combining operators
 * Forking operators
 * Maths operators
 * Other

In this episode you will see example of different types of operators in action.


## Other

### view

The `view` operator prints the items emitted by a channel to the console appending a *new line* character to each item in the channel. For example:
~~~
Channel
      .of('1', '2', '3')
      .view()
~~~
{: .source}

It prints:

~~~
1
2
3
~~~
{: .source}

An optional *closure* `{}` parameter can be specified to customise how items are printed. For example:

#### Closures

Briefly, a closure is a block of code that can be passed as an argument to a function. Thus, you can define a chunk of code and then pass it around as if it were a string or an integer.

~~~
Channel.of('1', '2', '3').view({ "chr$it" })
~~~
{: .source}

It prints:

~~~
chr1
chr2
chr3
~~~
{: .output}

## Filtering operators

We can reduce the number of items in a channel by using filtering operators.

Here we will use the `filter` operator on the chr_ch channel specifying the type qualifier `Number` so that only numbers are returned. We will then use the `view` opertor to inspect the contents.

~~~
chr_ch = Channel.of( 1..21, 'X', 'Y' )
autosomes_ch =chr_ch.filter( Number )
autosomes_ch.view()
~~~
{: .source}


Operators can also be chained together.

The previous example could be written like.

~~~
chr_ch = Channel.of( 1..21, 'X', 'Y' )
            .filter( Number )
            .view()
~~~

~~~
chr_ch = Channel.of( 1..21, 'X', 'Y' )
autosomes_ch =chr_ch.filter( Number ).filter({ it % 2 == 0 }).view()
~~~
{: .source}

> # Closures
> In the above example the filter condition is wrapped in curly  brackets, instead of round brackets, since it specifies a closure as the operator’s argument. This just is a language syntax-sugar for filter({ it % 2 == 0 } )
{: .callout}

## Transforming operators

As the name suggests transforming operators are used to transform the items emitted by a channel to new values.

### map

The `map` operator applies a function of your choosing to every item emitted by a channel, and returns the items so obtained as a new channel. The function applied is called the mapping function and is expressed with a closure as shown in the example below:

~~~
Channel.of( 'chr1', 'chr2' )
    .map ({ it.replaceAll("chr","") })
    .view()
~~~
{: .source}

Here the map function uses the string function `replaceAll` to remove the chr prefix from each element.

~~~
1
2
~~~
{: .output}

A map can associate to each element a tuple containing any data as needed.

~~~
Channel
    .fromPath( 'data/ggal/*.fq' )
    .map ({ file -> [file, file.countFastq()] })
    .view ({ file, numreads -> "file $file contains $numreads reads" })
~~~
{: .source}



~~~

~~~
{: .output}

> ## map operator
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



> ## channel names separator
> > Note the use in this example of curly brackets and the `;` as channel names separator. This is needed because the actual parameter of into is a closure which defines the target channels to which the source one is connected.

###  flatten

The `flatten` operator transforms a channel in such a way that every item in a list or tuple is flattened so that each single entry is emitted as a sole element by the resulting channel.

~~~
list1 = [1,2,3]
list2 = [4, 5, 6]

Channel
    .of(list1, list2)
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

The `collect` operator collects all the items emitted by a channel to a list and return the resulting object as a sole emission. This can be extremely useful when combing the results from output of a process.

~~~
Channel
    .of( 1, 2, 3, 4 )
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
     .of( [1,'A'], [1,'B'], [2,'C'], [3, 'B'], [1,'C'], [2, 'A'], [3, 'D'] )
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


> ## Group Tuple
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


## Combing Operators

Combing operators allow you to merge channels.

### mix

The `mix` operator combines the items emitted by two (or more) channels into a single channel.
~~~
ch1 = Channel.of( 1,2,3 )
ch2 = Channel.of( 'X','Y' )
ch3 = Channel.of( 'mt' )

ch1 .mix(ch2,ch3).view()
~~~
{: .source}

~~~
1
2
3
X
Y
mt
~~~
{: .output}

The items in the resulting channel have the same order as in respective original channel, however there’s no guarantee that the element of the second channel are append after the elements of the first.



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
{: .output}

## Forking operators

Forking operators.

### into

The `into` operator connects a source channel to two or more target channels in such a way the values emitted by the source channel are copied to the target channels. For example:


~~~
Channel
     .of( 'chr1', 'chr2', 'chr3' )
     .into{ ch1; ch2 }

ch1.view{ "ch1 emits: " + it }
ch2.view{ "ch2 emits: " + it }
~~~
{: .source}

Produces.

~~~
ch1 emits: chr1
ch1 emits: chr2
ch2 emits: chr1
ch1 emits: chr3
ch2 emits: chr2
ch2 emits: chr3
~~~
{: .output}

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

Produces.

~~~
N E X T F L O W  ~  version 20.10.0
Launching `hello.nf` [prickly_swanson] - revision: e33ef1057f
1 is small
2 is small
3 is small
40 is large
50 is large
~~~
{: .output}

> ## multi-channel object
> The branch operator returns a multi-channel object i.e. a variable that holds more than one channel object.
{: .callout}

## Maths operators

The maths operators allows you to apply simple math function  on channels.

### count

The count operator creates a channel that emits a single item: a number that represents the total number of items emitted by the source channel. For example:

~~~
Channel
    .of(1..21,'X','Y')
    .count()
    .view()
~~~

## Splitting operators

These operators are used to split items emitted by channels into chunks that can be processed by downstream operators or processes.

The available splitting operators are:

* [splitCsv](https://www.nextflow.io/docs/latest/operator.html#splitcsv)

* [splitFasta](https://www.nextflow.io/docs/latest/operator.html#splitfasta)

* [splitFastq](https://www.nextflow.io/docs/latest/operator.html#splitfastq)

* [splitText](https://www.nextflow.io/docs/latest/operator.html#splittext)

### splitCsv

The splitCsv operator allows you to parse text items emitted by a channel, that are formatted using the CSV format, and split them into records or group them into list of records with a specified length. This is useful when you want to create a samplesheet.

In the simplest case just apply the splitCsv operator to a channel emitting a CSV formatted text files or text entries. For example:

~~~
csv_ch=Channel
    .of('sample_id,fastq_1,fastq_2\ngut1,data/ggal/gut_1.fq,gut_2.fq\nliver_1,data/ggal/liver_1.fq,liver_2.fq')
    .splitCsv()
csv_ch.view()
~~~
{: .source}

The above example shows hows CSV text is parsed and is split into single rows. Values can be accessed by its column index in the row object.

~~~
csv_ch.view({it[0]})
~~~
{: .source}

When the CSV begins with a header line defining the column names, you can specify the parameter `header: true` which allows you to reference each value by its name, as shown in the following example:

> ## splitCsv
>
> 1. Modify the above script to print the second row.
> 2. Modify the above script to include `header: true` and view all data
>
> > ## Solution
> > ~~~~
> > 1. csv_ch.view({it[0]})
> > 2. .splitCsv(header:true)
> > 2. csv_ch.view()
> > ~~~

> > ~~~
> {: .solution}
{: .challenge}

## More resources

Check the operators [documentation](https://www.nextflow.io/docs/latest/operator.html) on Nextflow web site.

{: .output}
{% include links.md %}
