---
title: "Operators"
teaching: 30
exercises: 10
questions:
- "How can I change the contents of a channel?"
- "How do I perform operations such as filtering on channels?"
- "What are the different kinds of operators?"
- "How do I combine operators?"
- "How do I process a CSV file?"
objectives:
- "Describe what Nextflow operators are."
- "Modify the contents/elements of a channel using an operator"
- "Perform filtering and combining operations on a channel. "
- "Use the `splitCsv` operator to parse text items emitted by a channel, that are formatted using the CSV format, ."

keypoints:
- "Nextflow *operators* are methods that allow you to modify, set or view channels."
- "Operators can be separated in to several groups; filtering , transforming , splitting , combining , forking and Maths operators"
- "To use an opertaor use the dor notation after the Channel object e.g. `my_ch.view()`."
- "You can parse text items emitted by a channel, that are formatted using the CSV format,  using the `splitCsv` operator."
---

# Operators

In the previous Channels episode we learnt how to create Nextflow channels to enable us to pass data and values around our workflow. If we want to modify the contents or behaviour of a channel Nextflow provides methods called `Operators`. We have previously used the `view` operator to view the contents of a channel. There are many more operator methods that can be applied to Nextflow channels that can be usefully separated into several groups:


 * Filtering operators
 * Transforming operators
 * Splitting operators
 * Combining operators
 * Forking operators
 * Maths operators
 * Other

In this episode you will see examples, and get to use different types of operators.

# Using Operators

To use an operator add , use a dot `.` , followed by the operator name and brackets `()` on a channel object.

~~~
channel_obj.<operator>()
~~~
{: .language-groovy }

### view

The `view` operator prints the items emitted by a channel to the console appending a *new line* character to each item in the channel.

~~~
ch= channel.of('1', '2', '3')
ch.view()
~~~
{: .language-groovy }


To make code more readable we can spit the operators over several lines.
~~~
channel
      .of('1', '2', '3')
      .view()
~~~
{: .language-groovy }

prints:

~~~
1
2
3
~~~
{: .language-groovy }


#### Closures

An optional *closure* `{}` parameter can be specified to customise how items are printed.

Briefly, a closure is a block of code that can be passed as an argument to a function. In this way you can define a chunk of code and then pass it around as if it were a string or an integer. By default the parameters for a closure are specified with the groovy keyword `$it`.

For example here we use the the `view` operator and apply a closure to it, to add a chr prefix to each element of the channel using string interpolation.

~~~
channel
  .of('1', '2', '3')
  .view({ "chr$it" })
~~~
{: .language-groovy }

It prints:

~~~
chr1
chr2
chr3
~~~
{: .output}

## Filtering operators

We can reduce the number of items in a channel by using filtering operators.

The `filter` operator allows you to get only the items emitted by a channel that satisfy a condition and discarding all the others. The filtering condition can be specified by using either a

* regular expression,
* a literal value,
* a type qualifier (i.e. a Java class), e.g. Number, String, Boolean
* or any boolean statement.

#### type qualifier

Here we will use the `filter` operator on the `chr_ch` channel specifying the  type qualifier `Number` so that only numeric items are returned. We will then use the `view` operator to print the contents.

~~~
chr_ch = channel.of( 1..22, 'X', 'Y' )
autosomes_ch =chr_ch.filter( Number )
autosomes_ch.view()
~~~
{: .language-groovy }

We can chained together multiple operators using a `.` .

The previous example could be rewritten like.

~~~
chr_ch = channel
  .of( 1..22, 'X', 'Y' )
  .filter( Number )
  .view()
~~~
{: .language-groovy }

#### regular expression

To filter by a regular expression you have to do is to put `~` right in front of the string literal regular expression (e.g. `~"(^[Nn]extflow)"` or using slashy strings. `~/^[Nn]extflow`).

The following example shows how to filter a channel by using a regular expression `~/^1.*/` that returns only strings that begin with 1:

~~~
chr_ch = channel
  .of( 1..22, 'X', 'Y' )
  .filter(~/^1.*/)
  .view()
~~~
{: .language-groovy }

~~~
1
10
11
12
13
14
15
16
17
18
19
~~~
{: .output}

#### Boolean statement

A filtering condition can be defined by using any a Boolean expression,. A Boolean expression, is expressed by a closure,`{}`, returning a boolean value. For example the following fragment shows how to combined a filter for a type qualifier `Number` and then combining another filter operator using a Boolean expression to emitting numbers less the 5:

~~~
channel
  .of( 1..22, 'X', 'Y' )
  .filter(Number)
  .filter {it<5}
  .view()
~~~
{: .language-groovy }

~~~
1
2
3
4
~~~
{: .output }

> # Closures
> In the above example the filter condition is wrapped in curly brackets, instead of round brackets, since it specifies a closure as the operator’s argument. This just is a language syntax-sugar for filter({ it<5})
{: .callout}

####  literal value

Finally if we only want to include elements of a specific value we can specify a literal value. In the example below we use the literal value `X` to filter the channel for only those elements containing the value `X`.

~~~
channel
  .of( 1..22, 'X', 'Y' )
  .filter('X')
  .view()
~~~

~~~
X
~~~
{: .output }


> ## Filter a channel
>
> Add the boolean statement filter `filter({ it % 2 == 0 })` to the Nextflow script below to view only the even numbered chromosomes.
> ~~~
> chr_ch = channel
>  .of( 1..22, 'X', 'Y' )
>  .view()
> ~~~
> {: .language-groovy }
> > ## Solution
> >
> > ~~~
> > chr_ch = channel
> >   .of( 1..22, 'X', 'Y' )
> >   .filter( Number )
> >   .filter({ it % 2 == 0 })
> >   .view()
> > ~~~    
> > {: .language-groovy }
> {: .solution}
{: .challenge}

## Transforming operators

If we we want to modify the items emitted by a channel we use transforming operators.

### map

The `map` operator applies a function of your choosing to every item emitted by a channel, and returns the items so obtained as a new channel. The function applied is called the mapping function and is expressed with a closure as shown in the example below:

~~~
channel.of( 'chr1', 'chr2' )
    .map ({ it.replaceAll("chr","") })
    .view()
~~~
{: .language-groovy }

Here the map function uses the groovy string function `replaceAll` to remove the chr prefix from each element.

~~~
1
2
~~~
{: .output}

We can also use the `map` operator to associate a tuple to each element.

In the example below we use the `map` operator to transform a channel containing fastq files to a new channel containing a tuple with the fastq file and the number of reads in the fastq file. We use the `countFastq` file method to count the number of records in a FASTQ formatted file.

We can change the default name of the closure parameter keyword from `it` to a more meaningful name `file` using  `->`. When we have multiple parameters we can specify the keywords at the start of the closure, e.g. `file, name ->`.

~~~
channel
    .fromPath( 'data/yeast/reads/*.fq.gz' )
    .map ({ file -> [file, file.countFastq()] })
    .view ({ file, numreads -> "file $file contains $numreads reads" })
~~~
{: .language-groovy }

This would produce.

~~~
file data/yeast/reads/ref1_2.fq.gz contains 14677 reads
file data/yeast/reads/etoh60_3_2.fq.gz contains 26254 reads
file data/yeast/reads/temp33_1_2.fq.gz contains 20593 reads
file data/yeast/reads/temp33_2_1.fq.gz contains 15779 reads
file data/yeast/reads/ref2_1.fq.gz contains 20430 reads
[..truncated..]
~~~
{: .output}

> ## map operator
>
> Add a map operator to the Nextflow script below to transform the contents into a tuple with the file's basename, using the `.baseName`, method. Finally print the resulting channel.
> ~~~
>  channel
>  .fromPath( 'data/yeast/reads/*.fq.gz' )
>  .view()
> ~~~
> {: .language-groovy }
> > ## Solution
> >
> > ~~~
> > channel
> >   .fromPath( 'data/yeast/reads/*.fq.gz' )
> >   .map ({file -> [ file.name, file.baseName ]})
> >   .view({name, file -> "file's basename: $name"})
> > ~~~
> > {: .language-groovy }
> {: .solution}
{: .challenge}


###  flatten

The `flatten` operator transforms a channel in such a way that every item in a `list` or `tuple` is flattened so that each single entry is emitted as a sole element by the resulting channel.

~~~
list1 = [1,2,3]

println("without flatten:")
channel
    .of(list1)
    .view()

println("with flatten:")
channel
    .of(list1)
    .flatten()
    .view()

~~~
{: .language-groovy }

The above snippet prints:

~~~
without flatten:
[1, 2, 3]
with flatten:
1
2
3
~~~
{: .output}

### collect

The reverse of the `flatten` operator is `collect`. The `collect` operator collects all the items emitted by a channel to a list and return the resulting object as a sole emission. This can be extremely useful when combing the results from the output of multiple processes, or a single process run multiple times.

~~~
channel
    .of( 1, 2, 3, 4 )
    .collect()
    .view()
~~~
{: .language-groovy }


It prints a single value:

~~~
[1,2,3,4]
~~~
{: .output}

The result of the collect operator is a `value channel` and can be used multiple times.

### groupTuple

The `groupTuple` operator collects `tuples` or `lists` of values by grouping together the channel elements that share the same key. Finally it emits a new tuple object for each distinct key collected.

For example.

~~~
channel
     .of( ['wt','wt_1.fq'], ['wt','wt_1.fq'], ["mut",'mut_1.fq'], ['mut', 'mut_2.fq'] )
     .groupTuple()
     .view()
~~~~     
{: .language-groovy }

~~~
[wt, [wt_1.fq, wt_1.fq]]
[mut, [mut_1.fq, mut_2.fq]]
~~~
{: .output }

If we know the number of items to be grouped we can use the `groupTuple` size parameter.
When the specified size is reached, the tuple is emitted. By default incomplete tuples (i.e. with less than size grouped items) are discarded (default).

For example.

~~~
channel
     .of( ['wt','wt_1.fq'], ['wt','wt_1.fq'], ["mut",'mut_1.fq'])
     .groupTuple(size:2)
     .view()
~~~~     
{: .language-groovy }

outputs,

~~~
[wt, [wt_1.fq, wt_1.fq]]
~~~
{: .output}

This operator is useful to process altogether all elements for which there’s a common property or a grouping key.

> ## Group Tuple
>  ~~~
>  channel.fromPath('data/yeast/reads/*.fq.gz')
>         .view()
> ~~~
> {: .language-groovy }
> Modify the Nextflow script above to use the `map` operator to associate to each file the name prefix using  the closure.
> ~~~
> `{ file -> [ file.name.split('_')[0], file ] }`
> ~~~
> {: .language-groovy }
> Finally group together all files having the same common prefix using the `groupTuple` operator.
> > ## Solution
> >
> > ~~~
> > channel.fromPath('data/yeast/reads/*.fq.gz')
> >     .map { file -> [ file.name.split('_')[0], file ] }
> >     .groupTuple()
> >     .view()
> > ~~~
> > {: .language-groovy }
> {: .solution}
{: .challenge}


## Combining Operators

Combining operators allow you to merge channels together. This can be useful when you want to combine the output channels from multiple processes to perform another task such as joint QC.

### mix

The `mix` operator combines the items emitted by two (or more) channels into a single channel.
~~~
ch1 = channel.of( 1,2,3 )
ch2 = channel.of( 'X','Y' )
ch3 = channel.of( 'mt' )

ch1 .mix(ch2,ch3).view()
~~~
{: .language-groovy }

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
reads1_ch = channel.from(['wt', 'wt_1.fq'], ['mut','mut_1.fq'])
reads2_ch= channel.from(['wt', 'wt_2.fq'], ['mut','mut_2.fq'])
reads1_ch.join(reads2_ch).view()
~~~
{: .language-groovy }

The resulting channel emits:
~~~
[wt, wt_1.fq, wt_2.fq]
[mut, mut_1.fq, mut_2.fq]
~~~
{: .output}

## Forking operators

Forking operators split a single channel into multiple channels.

### into

The `into` operator connects a source channel to two or more target channels in such a way the values emitted by the source channel are copied to the target channels. For example:


~~~
channel
     .of( 'chr1', 'chr2', 'chr3' )
     .into{ ch1; ch2 }

ch1.view()
ch2.view()
~~~
{: .language-groovy }

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


## Maths operators

The maths operators allows you to apply simple math function  on channels.

### count

The `count` operator creates a channel that emits a single item: a number that represents the total number of items emitted by the source channel. For example:

~~~
channel
    .of(1..22,'X','Y')
    .count()
    .view()
~~~
{: .language-groovy }

~~~
24
~~~
{: .output }

## Splitting operators

These operators are used to split items emitted by channels into chunks that can be processed by downstream operators or processes.

The available splitting operators are:

|[splitCsv](https://www.nextflow.io/docs/latest/operator.html#splitcsv)|The splitCsv operator allows you to parse text items emitted by a channel, that are formatted using the CSV format, and split them into records or group them into list of records with a specified length.|
|[splitFasta](https://www.nextflow.io/docs/latest/operator.html#splitfasta)|The splitFasta operator allows you to split the entries emitted by a channel, that are formatted using the FASTA format. It returns a channel which emits text item for each sequence in the received FASTA content.|
|[splitFastq](https://www.nextflow.io/docs/latest/operator.html#splitfastq)|The splitFastq operator allows you to split the entries emitted by a channel, that are formatted using the FASTQ format. It returns a channel which emits a text chunk for each sequence in the received item.|
|[splitText](https://www.nextflow.io/docs/latest/operator.html#splittext)|The splitText operator allows you to split multi-line strings or text file items, emitted by a source channel into chunks containing n lines, which will be emitted by the resulting channel.|

### splitCsv

The `splitCsv` operator allows you to parse text items emitted by a channel, that are formatted using the CSV format, and split them into records or group them into list of records with a specified length. This is useful when you want to use a sample sheet.

In the simplest case just apply the `splitCsv` operator to a channel emitting a CSV formatted text files or text entries. For example:

For the CSV file `samples.csv`.

~~~
cat data/yeast/samples.csv
~~~
{: .language-bash }


~~~
sample_id,fastq_1,fastq_2
ref1,data/yeast/reads/ref1_1.fq.gz,data/yeast/reads/ref1_2.fq.gz
ref2,data/yeast/reads/ref2_1.fq.gz,data/yeast/reads/ref2_2.fq.gz
~~~
{: .output }

We can use the `splitCsv()` operator to split the channel contaning a CSV file into three elements.

~~~
csv_ch=channel
    .fromPath('data/yeast/samples.csv')
    .splitCsv()
csv_ch.view()
~~~
{: .language-groovy }

~~~
[sample_id, fastq_1, fastq_2]
[ref1, data/yeast/reads/ref1_1.fq.gz, data/yeast/reads/ref1_2.fq.gz]
[ref2, data/yeast/reads/ref2_1.fq.gz, data/yeast/reads/ref2_2.fq.gz]
~~~
{: .output }

The above example shows hows the CSV file `samples.csv` is parsed and is split into three elements.

#### Accessing values

Values can be accessed by its positional index using the square brackets syntax`[index]`. So to access the first column you would use `[0]` as shown in the following example:

~~~
csv_ch=channel
    .fromPath('data/yeast/samples.csv')
    .splitCsv()
csv_ch
  .view({it[0]})
~~~
{: .language-groovy }

~~~
sample_id
ref1
ref2
~~~
{: .output }


#### Column headers

When the CSV begins with a header line defining the column names, you can specify the parameter `header: true` which allows you to reference each value by its name, as shown in the following example:

~~~
csv_ch=channel
    .fromPath('data/yeast/samples.csv')
    .splitCsv(header:true)
csv_ch.view({it.fastq_1})
~~~
{: .language-groovy }

~~~
data/yeast/reads/ref1_1.fq.gz
data/yeast/reads/ref2_1.fq.gz
~~~
{: .output}


> ## Parse a CSV file
>
>  Modify the Nextflow script to print the first column `sample_id`.
>  ~~~
> csv_ch=channel
>    .fromPath('data/yeast/samples.csv')
>    .splitCsv(header:true)
> csv_ch.view({it.fastq_1})
>  ~~~
> {: .language-groovy }
> > ## Solution
> > ~~~~
> >  csv_ch=channel
> >         .fromPath('data/yeast/samples.csv')
> >         .splitCsv(header:true)
> >
> > csv_ch.view({it.sample_id})
> > ~~~
> > {: .language-groovy }
> {: .solution}
{: .challenge}

## More resources

Theses are just a few of the operators see the operators [documentation](https://www.nextflow.io/docs/latest/operator.html) on Nextflow web site for more deatils.

{: .output}
{% include links.md %}
