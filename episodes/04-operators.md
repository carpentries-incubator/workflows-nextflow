---
title: "Operators"
teaching: 0
exercises: 0
questions:
- "Key question (FIXME)"
objectives:
- "First learning objective. (FIXME)"
keypoints:
- "First key point. Brief Answer to questions. (FIXME)"
---

# Operators

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
~~`
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

> ## Challenge Title
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

{% include links.md %}

