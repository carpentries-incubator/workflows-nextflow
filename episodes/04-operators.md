---
title: "Operators"
teaching: 30
exercises: 10
questions:
- "What are Nextflow operators?"
- "How do you perform operations such as filtering on channels?"
- "What are the different kinds of operators?"
- "How do you combine operators?"
- "How dcan I used opertors to process a .csv file?"
objectives:
- "Describe what Nextflow operators are."
- "Modify the contents/elements of a channel using a operator"
- "Perform filtering and combining operations on a channel. "
- "Use the `splitCsv` operator to read in a .csv file."

keypoints:
- "Nextflow *operators* are methods that allow you to modify, set or view channels."
- "Operators can be separated in to  groups; filtering , transforming , splitting , combining , forking and Maths operators"
- "You can split csv file using the splitCsv operator."
---

# Operators

 *Operators* are methods that can be applied to channel objects. We have previously used the `view` operator to view the contents of a channel. There are many more operator methods that can be applied to channel that can be separated in to several groups:


 * Filtering operators
 * Transforming operators
 * Splitting operators
 * Combining operators
 * Forking operators
 * Maths operators
 * Other

In this episode you will see example, and get to use different types of operators.

# Using Operators

To use an operator add the operator name to channel object.

~~~
ch= Channel.of('1', '2', '3')
ch.view()
~~~
{: .language-groovy }

### view

The `view` operator prints the items emitted by a channel to the console appending a *new line* character to each item in the channel.

To make code more readable we can spit the operators over several lines.

~~~
Channel
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
Briefly, a closure is a block of code that can be passed as an argument to a function. In this way you can define a chunk of code and then pass it around as if it were a string or an integer. By default the arguments for a closure are specified with the variable `$it`.

~~~
Channel
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

Here we will use the `filter` operator on the `chr_ch` channel specifying the argument `Number` so that only numeric items are returned. We will then use the `view` operator to inspect the contents.

~~~
chr_ch = Channel.of( 1..21, 'X', 'Y' )
autosomes_ch =chr_ch.filter( Number )
autosomes_ch.view()
~~~
{: .language-groovy }


Operators can  be chained together.

The previous example could be written like.

~~~
chr_ch = Channel
  .of( 1..21, 'X', 'Y' )
  .filter( Number )
  .view()
~~~

Here we are using a closure as an argument to the `filter` method.

~~~
chr_ch = Channel
  .of( 1..21, 'X', 'Y' )
  .filter( Number )
  .filter({ it % 2 == 0 })
  .view()
~~~
{: .language-groovy }

> # Closures
> In the above example the filter condition is wrapped in curly  brackets, instead of round brackets, since it specifies a closure as the operator’s argument. This just is a language syntax-sugar for filter({ it % 2 == 0 } )
{: .callout}

## Transforming operators

Transforming operators are used to modify the items emitted by a channel to new values.

### map

The `map` operator applies a function of your choosing to every item emitted by a channel, and returns the items so obtained as a new channel. The function applied is called the mapping function and is expressed with a closure as shown in the example below:

~~~
Channel.of( 'chr1', 'chr2' )
    .map ({ it.replaceAll("chr","") })
    .view()
~~~
{: .language-groovy }

Here the map function uses the string function `replaceAll` to remove the chr prefix from each element.

~~~
1
2
~~~
{: .output}

The map method can associate to each element a tuple containing any data as needed.

~~~
Channel
    .fromPath( 'data/reads/*.fq.gz' )
    .map ({ file -> [file, file.countFastq()] })
    .view ({ file, numreads -> "file $file contains $numreads reads" })
~~~
{: .language-groovy }



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
> Use `fromPath` to create a channel emitting the fastq files matching the pattern `data/yeast/reads/*.fq.gz`, then chain with a map to return a pair containing the file name and basename using the `.baseName` method. Finally print the resulting channel.
>
> > ## Solution
> >
> > ~~~
> > Channel
> >   .fromPath( 'data/ggal/*.fq' )
> >   .map ({file -> [ file.name, file.baseName ]})
> >   .view({name, file -> "> file: $name"})
> > ~~~    
> {: .solution}
{: .challenge}


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
{: .language-groovy }
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

The `collect` operator collects all the items emitted by a channel to a list and return the resulting object as a sole emission. This can be extremely useful when combing the results from the output of multiple processes, or a single process run multiple times.

~~~
Channel
    .of( 1, 2, 3, 4 )
    .collect()
    .view()
~~~
{: .language-groovy }


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
     .of( ['wt','wt_1.fq'], ['wt','wt_1.fq'], ["mut",'mut_1.fq'], ['mut', 'mut_2.fq'] )
     .groupTuple()
     .view()
~~~~     
{: .language-groovy }

It shows:
~~~
[wt, [wt_1.fq, wt_1.fq]]
[mut, [mut_1.fq, mut_2.fq]]
~~~

This operator is useful to process altogether all elements for which there’s a common property or a grouping key.

Exercise


> ## Group Tuple
>
> 1. Use `fromPath` to create a channel emitting the fastq files matching the pattern `data/yeast/*.fq.gz`,
> 2. Use a map to associate to each file the name prefix.
> 3. Finally group together all files having the same common prefix.
>
> > ## Solution
> >
> > ~~~
> > Channel.fromPath('data/yeast/reads/*.fq.gz')
> >     .map { file -> [ file.name.split('_')[0], file ] }
> >     .groupTuple()
> >     .view()
> > ~~~
> {: .solution}
{: .challenge}


## Combining Operators

Combining operators allow you to merge channels together. This can be useful when you want to combine the output channels from multiple processes to perform another task such as joint QC.

### mix

The `mix` operator combines the items emitted by two (or more) channels into a single channel.
~~~
ch1 = Channel.of( 1,2,3 )
ch2 = Channel.of( 'X','Y' )
ch3 = Channel.of( 'mt' )

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
reads1_ch = Channel.from(['wt', 'wt_1.fq'], ['mut','mut_1.fq'])
reads2_ch= Channel.from(['wt', 'wt_2.fq'], ['mut','mut_2.fq'])
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

Forking operators.

### into

The `into` operator connects a source channel to two or more target channels in such a way the values emitted by the source channel are copied to the target channels. For example:


~~~
Channel
     .of( 'chr1', 'chr2', 'chr3' )
     .set{ ch1; ch2 }

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
{: .language-groovy }

The above example shows hows CSV text is parsed and is split into single rows. Values can be accessed by its column index in the row object.

~~~
csv_ch.view({it[0]})
~~~
{: .language-groovy }

When the CSV begins with a header line defining the column names, you can specify the parameter `header: true` which allows you to reference each value by its name, as shown in the following example:

~~~
csv_ch=Channel
    .of('sample_id,fastq_1,fastq_2\ngut1,data/ggal/gut_1.fq,gut_2.fq\nliver_1,data/ggal/liver_1.fq,liver_2.fq')
    .splitCsv(header:true)
csv_ch.view({it.fastq_1})
~~~

> ## splitCsv
>
> 1. Modify the above script to print the first column `sample_id`.
>
> > ## Solution
> > ~~~~
> >  csv_ch.splitCsv(header:true).view({it$sample_id})
> > ~~~

> > ~~~
> {: .solution}
{: .challenge}

## More resources

Theses are just a few of the operators see the operators [documentation](https://www.nextflow.io/docs/latest/operator.html) on Nextflow web site for more deatils.

{: .output}
{% include links.md %}
