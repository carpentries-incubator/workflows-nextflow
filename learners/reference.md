---
title: 'Groovy syntax'
---

## Groovy syntax

Here are some fundamental
concepts of the Groovy language.

- Single line comments (lines not meant to be executed) start
  with `//`. Multi-line comments are nested between `/*` and `*/` tags.
  ```groovy
  // This is a single line comment. Everything after the // is ignored.
  /*
      Comments can also
      span multiple
      lines.
  */
  ```
- Variables are assigned using `=` and can have any value. Variables used
  inside a double quoted string are prefixed with a `$` to denote the
  variable should be interpolated (replace the variable with the variable's
  value). A variable is otherwise just used by it's
  name.
  ```groovy
  myvar = 1                                     // Integer
  myvar = -3.1499392                            // Floating point number
  myvar = false                                 // Boolean
  myvar = "Hello world"                         // String
  myvar = new java.util.Date()                  // Object - Abstract data structure
  message = "The file $file cannot be found!"   // Variable used inside a string.
  println message                               // Print variable
  ```
- Lists (also known as arrays) are defined using the square bracket `[]` notation.
  ```groovy
  emptylist = []                      // an empty list
  // Lists can contain duplicates, and the values can be of any type.
  mixedList = [1, 1, 'string', true, null, 5 as byte]
  mixedList.add("new value")          // adds "new value" to the end of mixedList
  println mixedList.size()            // prints the size
  println mixedList[0]                // prints 1
  ```
- Maps (also known as associative arrays) are defined using the `[:]` literal. They associate a unique string with a value, and are commonly referred to as key-value pairs.
  ```groovy
  emptyMap = [:]                      // an empty map
  mymap = [ name : "Steve", age: 43, likes: ['walks','cooking','coding']]
  mymap.age = 44                      // values can be accessed using map property notation.
  println mymap['name']               // alternatively values can be accessed using quoted keys.
  ```
- Groovy supports common control structures such as if/else tests,
  switch/case tests, for loops, and while loops.
  ```groovy
  // if / else test
  x = false
  if ( !x ) {
      return 1
  } else {
      return 0
  }
  
  // switch / case test
  switch (x) {
      case "found foo":
          result = "Got foo"
          // continues to following cases
  
      case Number:
          result = "Number"
          break          // stops at this case
  
      default:
          result = "default"
  }
  
  // for loop (other variants exist)
  String message = ''
  for (int i = 0; i < 5; i++) {
      message += 'Hi '
  }
  
  // while loop
  y = 5
  while ( y-- > 0 ){
      println "Are we there yet?"
  }
  println "We've arrived!"
  ```
- Closures are open, anonymous, blocks of code that can take arguments,
  return a value and be assigned to a variable. A closure definition
  follows the syntax `{ [closureParameters ->] statements }`. A closure
  is evaluated when it is called, returning the result of the last statement
  in the closure.
  ```groovy
  // Find all elements > 1
  above_one = [1, 2, 3].findAll { it -> it > 1 }
  
  // Closure in a string
  message = "The file ${ file.getName() } cannot be found!"
  ```

Groovy is very syntax-rich and supports many more operations. A full
description of Groovy semantics can be found in the [Groovy Documentation](https://groovy-lang.org/semantics.html).

## Glossary

| Term                     | Description                                                                                                                                                                                                                         | 
| ------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Cache                    | In Nextflow, caching refers to the mechanism of reusing the results of executed processes in subsequent runs. This feature reduces computation time by avoiding the re-execution of unchanged tasks.                                |
| Channel                  | Channels are the conduits through which data flows in a Nextflow workflow. They connect processes, allowing for the movement and transformation of data from one process to the next.                                               |
| Closure                  | In Nextflow, a closure is a block of code that can be assigned to a variable, passed as an argument, or executed.        |
| Dataflow programming     | Dataflow programming is a programming paradigm that models a program as a directed graph of the data flowing between operations.                                                                                                    |
| Directive                | Directives are annotations in a Nextflow script that provide execution instructions for processes. They control aspects like resource allocation, process retries, and container usage.                                             |
| Docker                   | Docker is an open-source platform used to automate the deployment of applications inside lightweight, portable containers that allow applications to run seamlessly in different computing environments.                            |
| Domain Specific Language | A domain-specific language (DSL) is a computer language specialized to a particular application domain. Nextflow is a language for computational workflows.                                                                         |
| DSL2                     | The second version of the Nextflow domain-specific language, introducing modular and reusable workflow components. DSL2 allows for more structured and maintainable pipeline development.                                           |
| Executor                 | The component responsible for managing the execution of processes in a Nextflow workflow. Executors handle task scheduling, execution, and resource allocation across different computing environments.                             |
| Groovy                   | The programming language used to write Nextflow scripts. Groovy is a powerful, optionally typed, and dynamic language, with static-typing and static compilation capabilities. It enhances Nextflow with scripting flexibility.     |
| Module                   | A reusable component in a Nextflow pipeline that encapsulates a specific function or set of functions. Modules can be combined and reused across different workflows, promoting modularity and code reusability.                    |
| Nextflow                 | A domain-specific language and toolkit for parallel and scalable computational workflows. Nextflow enables the scripting of complex data analysis pipelines in a reproducible and efficient manner.                                 |
| Nextflow workflow        | In Nextflow, a workflow is a composition of processes and dataflow logic (i.e., channels and operators).                                                                                                                            |
| Operators                | In Nextflow, operators are methods that allow the transformation, filtering, and combination of data streams. They facilitate the management and manipulation of data within a workflow.                                            |
| Parameters               | Variables defined in a Nextflow script or passed at runtime to customize the execution of a workflow. Parameters allow for the dynamic configuration of workflows based on input data or user-defined settings.                     |
| Pipeline                 | A sequence of processes orchestrated by Nextflow to perform a complex data analysis task. Pipelines are defined in Nextflow scripts and can be composed of multiple workflow blocks and processes.                                  |
| Process                  | A fundamental component in Nextflow representing a computational task. Each process defines a small part of the overall workflow, encapsulating the command-line tools, scripts, and computational resources required for execution.|
| Profile                  | A set of configuration options in Nextflow that can be activated to adapt the execution environment of a workflow. Profiles enable easy switching between different computing environments, like local, cloud, or cluster execution.|
| PublishDir               | A directive in a process block that specifies a directory where the output files of the process should be stored. It helps in organizing and accessing the results of workflow executions.                                          |
| Script Block             | A block within a Nextflow workflow that contains custom script code. Script blocks allow for the embedding of standard programming constructs and commands within the flow of a workflow.                                           |
| Tuple                    | An ordered, immutable list of elements. A tuple can be seen as an ordered collection of objects of different types. These objects do not necessarily relate to each other in any way, but collectively they will have some meaning. | 
| Workflow Block           | A segment of a Nextflow script defining a series of operations or processes to be executed.                                                                                                                                         |
| WorkDir                  | The working directory where Nextflow stores intermediate files generated by processes. It ensures reproducibility and efficient management of temporary data.                                                                       |



## External references

### Manuals

- [Nextflow documentation](https://www.nextflow.io/docs/latest/index.html)
- [nf-core](https://nf-co.re/)
- [nf-core pipelines](https://nf-co.re/pipelines)
- [A step-by-step usage guide to nf-core pipelines](https://nf-co.re/usage/introduction)
- [nf-core Troubleshooting guide](https://nf-co.re/usage/troubleshooting)
- [Example nf-core institutional config file](https://github.com/nf-core/configs/blob/master/conf/eddie.config)
- If you wanted to check if your institute has nf-core docs (for example a config file), they have a list of the systems [here](https://github.com/nf-coe/configs#documentation)
- Example of how you can usefully use `log.info` [here](https://github.com/nextflow-io/rnaseq-nf/blob/3b5b49f/main.nf#L41-L48_)
- [Nextflow have a slack space](https://www.nextflow.io/blog/2022/nextflow-is-moving-to-slack.html)
- [nf-core slack](https://nf-co.re/join)
- Here's some details on a nextflow issue relating to cancelling pipelines that are running. [Ability to stop running pipeline](https://github.com/nextflow-io/nextflow/issues/1441).
- [Groovy Documentation](https://groovy-lang.org/documentation.html)

### Papers

- [Using prototyping to choose a bioinformatics workflow management system](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008622)




