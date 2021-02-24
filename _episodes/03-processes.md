---
title: "Processes"
teaching: 0
exercises: 0
questions:
- "What is a Nextflow process?"
- "How do I create a Nextflow process?"
- "How do I input data into processes|"
- "How do I ouput data from a process?"
- "How do I sepcify conditions for a process in order for it to execute?"
- "What are process directives?"
- "How do i save output from a process?"
objectives:
- "Understand a Nextflow process."
- "Create a nextflow process."
- "Use values and files as inputs to a process."
- "Use the when declaration to define a condition for process execution."
- "Understand what process directives are." 
keypoints:
- "A Nextflow process is an independent task/step in a workflow"
- "Processes"
---


# Processes

A process is the basic Nextflow computing primitive to execute foreign function i.e. custom scripts or tools.

*Processes can be thought of as a particular task/steps in a workflow, e.g. an alignment step in RNA-Seq analysis. Processes and are independent of each other (don't require another processes to execute) and can not communicate/write to each other . It is the Channels that pass the data from each process to another, and we do this by having the processes define input and output which are Channels*

The process definition starts with keyword the process, followed by process name and finally the process body delimited by curly brackets. The process body must contain a string which represents the command or, more generally, a script that is executed by it.

A basic process looks like the following example:

~~~
process sayHello {
  """
  echo 'Hello world!'
  """
}
~~~
{: .source}

A process may contain five definition blocks, respectively: directives, inputs, outputs, when clause and finally the process script. The syntax is defined as follows:

~~~
process < name > {
  [ directives ]        
  input:                
  < process inputs >
  output:               
  < process outputs >
  when:                 
  < condition >
  [script|shell|exec]:  
  < user script to be executed >
}
~~~
{: .source}

* Zero, one or more process directives
* Zero, one or more process inputs
* Zero, one or more process outputs
* An optional boolean conditional to trigger the process execution
* The command to be executed


## Script

The script block is a string statement that defines the command that is executed by the process to carry out its task.

A process contains one and only one script block, and it must be the last statement when the process contains input and output declarations.

The script block can be a simple string or multi-line string. The latter simplifies the writing of non trivial scripts composed by multiple commands spanning over multiple lines. For example::


~~~
process example {
    script:
    """
    blastp -db /data/blast -query query.fa -outfmt 6 > blast_result
    cat blast_result | head -n 10 | cut -f 2 > top_hits
    blastdbcmd -db /data/blast -entry_batch top_hits > sequences
    """
}
~~~
{: .source}

By default the process command is interpreted as a **Bash** script. However any other scripting language can be used just simply starting the script with the corresponding [Shebang](https://en.wikipedia.org/wiki/Shebang_(Unix)) declaration. For example:

~~~
process pyStuff {
  script:
  """
  #!/usr/bin/env python

  x = 'Hello'
  y = 'world!'
  print "%s - %s" % (x,y)
  """
}
~~~

> This allows the compositing in the same workflow script of tasks using different programming languages which may better fit a particular job. However for large chunks of code is suggested to save them into separate files and invoke them from the process script.
{: .callout}


### Script parameters

Process script can be defined dynamically using variable values like in other string.

~~~
params.data = 'World'

process foo {
  script:
  """
  echo Hello $params.data
  """
}
~~~
{: .source}

> # String interpolation
> A process script can contain any string format supported by the Groovy programming language. This allows us to use string interpolation or multiline string as in the script above. Refer to [String interpolation](https://seqera.io/training/#_string_interpolation) for more information.
{: .callout}

> Since Nextflow uses the same Bash syntax for variable substitutions in strings, Bash environment variables need to be escaped using `\` character.
{: .callout}

~~~
process foo {
  script:
  """
  echo "The current directory is \$PWD"
  """
}
~~~
{: .source}

> ## Escape Bash
>
> Try to modify the above script using $PWD instead of \$PWD and check the difference.
.
>
> > ## Solution
> >
> > This is the body of the solution.
> {: .solution}
{: .challenge}

However this won’t allow any more the usage of Nextflow variables in the command script.

Another alternative is to use a `shell` statement instead of `script` which uses a different syntax for Nextflow variable: `!{..}`. This allow to use both Nextflow and Bash variables in the same script.

```
params.data = 'le monde'

process baz {
  shell:
  '''
  X='Bonjour'
  echo $X !{params.data}
  '''
}
```
{: .source}

### Conditional script

The process script can also be defined in a complete dynamic manner using a if statement or any other expression evaluating to string value. For example:

~~~
params.aligner = 'kallisto'

process foo {
  script:
  if( params.aligner == 'kallisto' )
    """
    kallisto --reads /some/data.fastq
    """
  else if( params.aligner == 'salmon' )
    """
    salmon --reads /some/data.fastq
    """
  else
    throw new IllegalArgumentException("Unknown aligner $params.aligner")
}
~~~
{: .source}

> ## Conditional Exercise
>
> Write a custom function that given the aligner name as parameter returns the command string to be executed. Then use this function as the process script body.
.
>
> > ## Solution
> >
> > This is the body of the solution.
> {: .solution}
{: .challenge}

## Inputs

Nextflow processes are isolated from each other but can communicate between themselves sending values through channels.

Inputs implicitly determine the dependency and the parallel execution of the process. The process execution is fired each time a new data is ready to be consumed from the input channel:

todo insert image channel-process

The input block defines which channels the process is expecting to receive inputs data from. You can only define one input block at a time and it must contain one or more inputs declarations.

The input block follows the syntax shown below:

~~~
input:
  <input qualifier> <input name> from <source channel>
~~~
{: .source}

### Input values

The val qualifier allows you to receive data of any type as input. It can be accessed in the process script by using the specified input name, as shown in the following example:

~~~
num = Channel.from( 1, 2, 3 )

process basicExample {
  input:
  val x from num

  """
  echo process job $x
  """
}
~~~
{: .source}

In the above example the process is executed three times, each time a value is received from the channel num and used to process the script. Thus, it results in an output similar to the one shown below:

~~~
process job 3
process job 1
process job 2
~~~
{: .output}

> ## Channel order
> The channel guarantees that items are delivered in the same order as they have been sent - but - since the process is executed in a parallel manner, there is no guarantee that they are processed in the same order as they are received.
{: .callout}

### Input files

The `file` qualifier allows the handling of file values in the process execution context. This means that Nextflow will stage it in the process execution directory, and it can be access in the script by using the name specified in the input declaration.

~~~
reads = Channel.fromPath( 'data/ggal/*.fq' )

process foo {
    input:
    file 'sample.fastq' from reads
    script:
    """
    your_command --reads sample.fastq
    """
}
~~~
{: .source}





{% include links.md %}

