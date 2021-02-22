---
title: "Channels"
teaching: 15
exercises: 0
questions:
- "Key question (FIXME)"
objectives:
- "First learning objective. (FIXME)"
keypoints:
- "First key point. Brief Answer to questions. (FIXME)"
---

# Channels

Channels are a key data structure of Nextflow that allows the implementation of reactive-functional oriented computational workflows based on the [Dataflow programming paradigm](https://en.wikipedia.org/wiki/Dataflow_programming).

They are used to logically connect tasks each other or to implement functional style data transformations.


## Channel types

Nextflow distinguish two different kinds of channels: queue channels and value channels.

### Queue channel

A queue channel is a asynchronous unidirectional FIFO queue which connects two processes or operators.

* What asynchronous means? That operations are non-blocking.

* What unidirectional means? That data flow from a producer to a consumer.

* What FIFO means? That the data is guaranteed to be delivered in the same order as it is produced.

A queue channel is implicitly created by process output definitions or using channel factories methods such as [Channel.from](https://www.nextflow.io/docs/latest/channel.html#from) or [Channel.fromPath](https://www.nextflow.io/docs/latest/channel.html#frompath).

Try the following snippets:


> ## View Channel contents
> ~~~
> ch = Channel.from(1,2,3)
> ~~~
> 	Use the built-in println function to print the ch variable.
>	  Apply the view method to the ch channel, therefore prints each item emitted by the channels.
>
> > ## Solution
> > ~~~
> > println(ch)  
> > ch.view() 
> > ~~~
> {: .solution}
{: .challenge}

> A queue channel can have one and exactly one producer and one and exactly one consumer.
{: .callout}





{% include links.md %}

