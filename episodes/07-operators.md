---
title: Operators
teaching: 30
exercises: 10
---

::::::::::::::::::::::::::::::::::::::: objectives

- Understand what Nextflow operators are.
- Modify the contents/elements of a channel using operators.
- Perform filtering and combining operations on a channel object.
- Use the `splitCsv` operator to parse the contents of  CSV file into a channel .

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::: questions

- How do I perform operations, such as filtering, on channels?
- What are the different kinds of operations I can perform on channels?
- How do I combine operations?
- How can I use a CSV file to process data into a Channel?

::::::::::::::::::::::::::::::::::::::::::::::::::

## Operators

In the Channels episode we learnt how to create Nextflow channels to enable us to pass data and values around our workflow. If we want to modify the contents or behaviour of a channel, Nextflow provides methods called `operators`. We have previously used the `view` operator to view the contents of a channel. There are many more operator methods that can be applied to Nextflow channels that can be usefully separated into several groups:

- **Filtering** operators: reduce the number of elements in a channel.
- **Transforming** operators: transform the value/data in a channel.
- **Splitting** operators: split items in a channel into smaller chunks.
- **Combining** operators: join channels together.
- **Maths** operators: apply simple math functions on channels.
- **Other**: such as the view operator.

In this episode you will see examples, and get to use different types of operators.

## Using Operators

To use an operator, the syntax is the channel name, followed by a dot `.` , followed by the operator name and brackets `()`.

```groovy 
channel_obj.<operator>()
```

#### view

The `view` operator prints the items emitted by a channel to the console appending a *new line* character to each item in the channel.

```groovy 
ch = channel.of('1', '2', '3')
ch.view()
```

We can also chain  together the channel factory method `.of` and the operator `.view()` using the
dot notation.

```groovy 
ch = channel.of('1', '2', '3').view()
```

To make code more readable we can split the operators over several lines.
The blank space between the operators is ignored and is solely for readability.

```groovy 
ch = channel
      .of('1', '2', '3')
      .view()
```

prints:

```groovy 
1
2
3
```

##### Closures

An optional *closure* `{}` parameter can be specified to customise how items are printed.

Briefly, a closure is a block of code that can be passed as an argument to a function. In this way you can define a chunk of code and then pass it around as if it were a string or an integer. By default the parameters for a closure are specified with the groovy keyword `$it` ('it' is for 'item').

For example here we use the the `view` operator and apply a closure to it, to add a `chr` prefix to each element of the channel using string interpolation.

```groovy 
ch = channel
  .of('1', '2', '3')
  .view({ "chr$it" })
```

It prints:

```output
chr1
chr2
chr3
```

**Note:** the `view()` operator doesn't change the contents of the channel object.

```groovy 
ch = channel
  .of('1', '2', '3')
  .view({ "chr$it" })

ch.view()  
```

```output
chr1
chr2
chr3
1
2
3
```

### Filtering operators

We can reduce the number of items in a channel by using filtering operators.

The `filter` operator allows you to get only the items emitted by a channel that satisfy a condition and discard all the others. The filtering condition can be specified by using either:

- a regular expression
- a literal value
- a data type qualifier, e.g. Number (any integer,float ...), String, Boolean
- or any boolean statement.

##### Data type qualifier

Here we use the `filter` operator on the `chr_ch` channel specifying the  data type qualifier `Number` so that only numeric items are returned. The Number data type includes both integers and floating point numbers.
We will then use the `view` operator to print the contents.

```groovy 
chr_ch = channel.of( 1..22, 'X', 'Y' )
autosomes_ch =chr_ch.filter( Number )
autosomes_ch.view()
```

To simplify the code we can chain multiple operators together, such as `filter` and `view` using a `.` .

The previous example could be rewritten like:
The blank space between the operators is ignored and is used for readability.

```groovy 
chr_ch = channel
  .of( 1..22, 'X', 'Y' )
  .filter( Number )
  .view()
```

```output 
1
2
3
4
5
6
7
8
9
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
20
21
22
```

##### Regular expression

To filter by a regular expression you have to do is to put `~` right in front of the string literal regular expression (e.g. `~"(^[Nn]extflow)"` or use slashy strings which replace the quotes with `/`. `~/^[Nn]extflow/`).

The following example shows how to filter a channel by using a regular expression `~/^1.*/` inside a slashy string, that returns only strings that begin with 1:

```groovy 
chr_ch = channel
  .of( 1..22, 'X', 'Y' )
  .filter(~/^1.*/)
  .view()
```

```output
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
```

##### Boolean statement

A filtering condition can be defined by using a Boolean expression described by a closure `{}` and returning a boolean value. For example the following fragment shows how to combine a filter for a type qualifier `Number` with another filter operator using a Boolean expression to emit numbers less than 5:

```groovy 
channel
  .of( 1..22, 'X', 'Y' )
  .filter(Number)
  .filter { it < 5 }
  .view()
```

```output 
1
2
3
4
```

:::::::::::::::::::::::::::::::::::::::::  callout

### Closures

In the above example we have removed the brackets around the filter condition e.g. `filter{ it<5}`, since it specifies a closure as the operator's argument. This is language short for `filter({ it<5})`


::::::::::::::::::::::::::::::::::::::::::::::::::

##### Literal value

Finally, if we only want to include elements of a specific value we can specify a literal value. In the example below we use the literal value `X` to filter the channel for only those elements containing the value `X`.

```
channel
  .of( 1..22, 'X', 'Y' )
  .filter('X')
  .view()
```

```output 
X
```

:::::::::::::::::::::::::::::::::::::::  challenge

### Filter a channel

Add two channel filters to the Nextflow script below to view only the even numbered chromosomes.

**Note:** The expression `it % 2`  produces the remainder of a division.

```groovy 
chr_ch = channel
 .of( 1..22, 'X', 'Y' )
 .view()
```

:::::::::::::::  solution

### Solution

```groovy 
chr_ch = channel
  .of( 1..22, 'X', 'Y' )
  .filter( Number )
  .filter({ it % 2 == 0 })
  .view()
```

```output 
2
4
6
8
10
12
14
16
18
20
22
```

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

### Modifying the contents of a channel

If we want to modify the items in a channel, we can use transforming operators.

#### map

#### Applying a function to items in a channel

The `map` operator applies a function of your choosing to every item in a channel, and returns the items so obtained as a new channel. The function applied is called the mapping function and is expressed with a closure `{}` as shown in the example below:

```groovy 
chr = channel
  .of( 'chr1', 'chr2' )
  .map ({ it.replaceAll("chr","") })

chr.view()
```

Here the map function uses the groovy string function `replaceAll` to remove the chr prefix from each element.

```output
1
2
```

We can also use the `map` operator to transform each element into a tuple.

In the example below we use the `map` operator to transform a channel containing fastq files to a new channel containing a tuple with the fastq file and the number of reads in the fastq file. We use the built in `countFastq` file method to count the number of records in a FASTQ formatted file.

We can change the default name of the closure parameter keyword from `it` to a more meaningful name `file` using  `->`. When we have multiple parameters we can specify the keywords at the start of the closure, e.g. `file, numreads ->`.

```groovy 
fq_ch = channel
    .fromPath( 'data/yeast/reads/*.fq.gz' )
    .map ({ file -> [file, file.countFastq()] })
    .view ({ file, numreads -> "file $file contains $numreads reads" })
```

This would produce.

```output
file data/yeast/reads/ref1_2.fq.gz contains 14677 reads
file data/yeast/reads/etoh60_3_2.fq.gz contains 26254 reads
file data/yeast/reads/temp33_1_2.fq.gz contains 20593 reads
file data/yeast/reads/temp33_2_1.fq.gz contains 15779 reads
file data/yeast/reads/ref2_1.fq.gz contains 20430 reads
[..truncated..]
```

We can then add a `filter` operator to only retain those fastq files with more than 25000 reads.

```groovy 
channel
    .fromPath( 'data/yeast/reads/*.fq.gz' )
    .map ({ file -> [file, file.countFastq()] })
    .filter({ file, numreads -> numreads > 25000})
    .view ({ file, numreads -> "file $file contains $numreads reads" })
```

```output
file data/yeast/reads/etoh60_3_2.fq.gz contains 26254 reads
file data/yeast/reads/etoh60_3_1.fq.gz contains 26254 reads
```

:::::::::::::::::::::::::::::::::::::::  challenge

### map operator

Add a `map` operator to the Nextflow script below to transform the contents into a tuple with the file and the file's name, using the `.getName` method. The `getName` method gives the filename. Finally `view` the channel contents.

```groovy 
 channel
 .fromPath( 'data/yeast/reads/*.fq.gz' )
 .view()
```

:::::::::::::::  solution

### Solution

```groovy 
ch = channel
  .fromPath( 'data/yeast/reads/*.fq.gz' )
  .map ({file -> [ file, file.getName() ]})
  .view({file, name -> "file's name: $name"})
```

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

#### Converting a list into multiple items

The `flatten` operator transforms a channel in such a way that every item in a `list` or `tuple` is flattened so that each single entry is emitted as a sole element by the resulting channel.

```groovy 
list1 = [1,2,3]
ch = channel
  .of(list1)
  .view()
```

```output
[1, 2, 3]
```

```groovy 
ch =channel
    .of(list1)
    .flatten()
    .view()

```

The above snippet prints:

```output
1
2
3
```

This is similar to the channel factory `Channel.fromList`.

#### Converting the contents of a channel to a single list item.

The reverse of the `flatten` operator is `collect`. The `collect` operator collects all the items emitted by a channel to a list and return the resulting object as a sole emission. This can be extremely useful when combining the results from the output of multiple processes, or a single process run multiple times.

```groovy 
ch = channel
    .of( 1, 2, 3, 4 )
    .collect()
    .view()
```

It prints a single value:

```output
[1,2,3,4]
```

The result of the collect operator is a `value` channel and can be used multiple times.

#### Grouping contents of a channel by a key.

The `groupTuple` operator collects `tuples` or `lists` of values by grouping together the channel elements that share the same key. Finally it emits a new tuple object for each distinct key collected.

For example.

```groovy 
ch = channel
     .of( ['wt','wt_1.fq'], ['wt','wt_2.fq'], ["mut",'mut_1.fq'], ['mut', 'mut_2.fq'] )
     .groupTuple()
     .view()
```

```output 
[wt, [wt_1.fq, wt_1.fq]]
[mut, [mut_1.fq, mut_2.fq]]
```

If we know the number of items to be grouped we can use the `groupTuple` `size` parameter.
When the specified size is reached, the tuple is emitted. By default incomplete tuples (i.e. with less than size grouped items) are discarded (default).

For example.

```groovy 
ch = channel
     .of( ['wt','wt_1.fq'], ['wt','wt_1.fq'], ["mut",'mut_1.fq'])
     .groupTuple(size:2)
     .view()
```

outputs,

```output
[wt, [wt_1.fq, wt_1.fq]]
```

This operator is useful to process altogether all elements for which there's a common property or a grouping key.

:::::::::::::::::::::::::::::::::::::::  challenge

### Group Tuple

```groovy 
channel.fromPath('data/yeast/reads/*.fq.gz')
       .view()
```

Modify the Nextflow script above to add the `map` operator to create a tuple with the name prefix as the key and the file as the value using the closure below.

```groovy 
{ file -> [ file.getName().split('_')[0], file ] }
```

Finally group together all files having the same common prefix using the `groupTuple` operator and `view` the contents of the channel.

:::::::::::::::  solution

### Solution

```groovy 
ch = channel.fromPath('data/yeast/reads/*.fq.gz')
    .map { file -> [ file.getName().split('_')[0], file ] }
    .groupTuple()
    .view()
```

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

### Merging Channels

Combining operators allows you to merge channels together. This can be useful when you want to combine the output channels from multiple processes to perform another task such as joint QC.

#### mix

The `mix` operator combines the items emitted by two (or more) channels into a single channel.

```groovy 
ch1 = channel.of( 1,2,3 )
ch2 = channel.of( 'X','Y' )
ch3 = channel.of( 'mt' )

ch4 = ch1.mix(ch2,ch3).view()
```

```output
1
2
3
X
Y
mt
```

The items emitted by the resulting mixed channel may appear in any order, regardless of which source channel they came from. Thus, the following example it could be a possible result of the above example as well.

```output
1
2
X
3
mt
Y
```

#### join

The `join` operator creates a channel that joins together the items emitted by two channels for which exists a matching key. The key is defined, by default, as the first element in each item emitted.

```groovy 
reads1_ch = channel
  .of(['wt', 'wt_1.fq'], ['mut','mut_1.fq'])
reads2_ch= channel
  .of(['wt', 'wt_2.fq'], ['mut','mut_2.fq'])
reads_ch = reads1_ch
  .join(reads2_ch)
  .view()
```

The resulting channel emits:

```output
[wt, wt_1.fq, wt_2.fq]
[mut, mut_1.fq, mut_2.fq]
```


### Maths operators

The maths operators allows you to apply simple math function on channels.

The maths operators are:

- count
- min
- max
- sum
- toInteger

#### Counting items in a channel

The `count` operator creates a channel that emits a single item: a number that represents the total number of items emitted by the source channel. For example:

```groovy 
ch = channel
    .of(1..22,'X','Y')
    .count()
    .view()
```

```output 
24
```

### Splitting items in a channel

Sometimes you want to split the content of a individual item in a channel, like a file or string, into smaller chunks that can be processed by downstream operators or processes e.g. items stored in a CSV file.

Nextflow has a number of splitting operators that can achieve this:

- [splitCsv](https://www.nextflow.io/docs/latest/operator.html#splitcsv): The splitCsv operator allows you to parse text items emitted by a channel, that are formatted using the CSV format, and split them into records or group them into list of records with a specified length.
- [splitFasta](https://www.nextflow.io/docs/latest/operator.html#splitfasta): The splitFasta operator allows you to split the entries emitted by a channel, that are formatted using the FASTA format. It returns a channel which emits a text item for each sequence in the received FASTA content.
- [splitFastq](https://www.nextflow.io/docs/latest/operator.html#splitfastq): The splitFastq operator allows you to split the entries emitted by a channel, that are formatted using the FASTQ format. It returns a channel which emits a text chunk for each sequence in the received item.
- [splitText](https://www.nextflow.io/docs/latest/operator.html#splittext): The splitText operator allows you to split multi-line strings or text file items, emitted by a source channel into chunks containing n lines, which will be emitted by the resulting channel.

#### splitCsv

The `splitCsv` operator allows you to parse text items emitted by a channel, that are formatted using the CSV format, and split them into records or group them into list of records with a specified length. This is useful when you want to use a sample sheet.

In the simplest case just apply the `splitCsv` operator to a channel emitting a CSV formatted text files or text entries. For example:

For the CSV file `samples.csv`.

```bash 
cat data/yeast/samples.csv
```

```output 
sample_id,fastq_1,fastq_2
ref1,data/yeast/reads/ref1_1.fq.gz,data/yeast/reads/ref1_2.fq.gz
ref2,data/yeast/reads/ref2_1.fq.gz,data/yeast/reads/ref2_2.fq.gz
```

We can use the `splitCsv()` operator to split the channel contaning a CSV file into three elements.

```groovy 
csv_ch=channel
    .fromPath('data/yeast/samples.csv')
    .splitCsv()
csv_ch.view()
```

```output 
[sample_id, fastq_1, fastq_2]
[ref1, data/yeast/reads/ref1_1.fq.gz, data/yeast/reads/ref1_2.fq.gz]
[ref2, data/yeast/reads/ref2_1.fq.gz, data/yeast/reads/ref2_2.fq.gz]
```

The above example shows hows the CSV file `samples.csv` is parsed and is split into three elements.

##### Accessing values

Values can be accessed by their positional indexes using the square brackets syntax`[index]`. So to access the first column you would use `[0]` as shown in the following example:

```groovy 
csv_ch=channel
    .fromPath('data/yeast/samples.csv')
    .splitCsv()
csv_ch
  .view({it[0]})
```

```output 
sample_id
ref1
ref2
```

##### Column headers

When the CSV begins with a header line defining the column names, you can specify the parameter `header: true` which allows you to reference each value by its name, as shown in the following example:

```groovy 
csv_ch=channel
    .fromPath('data/yeast/samples.csv')
    .splitCsv(header:true)
csv_ch.view({it.fastq_1})
```

```output
data/yeast/reads/ref1_1.fq.gz
data/yeast/reads/ref2_1.fq.gz
```

:::::::::::::::::::::::::::::::::::::::  challenge

### Parse a CSV file

Modify the Nextflow script to print the first column `sample_id`.

```groovy 
csv_ch=channel
  .fromPath('data/yeast/samples.csv')
```

:::::::::::::::  solution

### Solution

```groovy
 csv_ch=channel
        .fromPath('data/yeast/samples.csv')
        .splitCsv(header:true)

csv_ch.view({it.sample_id})
```

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

#### Tab delimited files

If you want to split a `tab` delimited file or file separated by another character use the `sep` parameter of the split `splitCsv` operator.

For examples,

```groovy 
Channel.of("val1\tval2\tval3\nval4\tval5\tval6\n")
  .splitCsv(sep: "\t")
  .view()
```

```output
[val1, val2, val3]
[val4, val5, val6]
```

### More resources

See the operators [documentation](https://www.nextflow.io/docs/latest/operator.html) on the Nextflow web site.



:::::::::::::::::::::::::::::::::::::::: keypoints

- Nextflow *operators* are methods that allow you to modify, set or view channels.
- Operators can be separated in to several groups; filtering , transforming , splitting , combining , forking and Maths operators
- To use an operator use the dot notation after the Channel object e.g. `my_ch.view()`.
- You can parse text items emitted by a channel, that are formatted using the CSV format,  using the `splitCsv` operator.

::::::::::::::::::::::::::::::::::::::::::::::::::


