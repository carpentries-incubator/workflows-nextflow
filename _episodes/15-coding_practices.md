---
title: "Nextflow coding practices"
teaching: 30
exercises: 15
questions:
- "How do I make my code readable?"
- "How do I make my code portable?"
- "How do I make my code maintainable?"
objectives:
- "Learn how to use whitespace and comments to improve code readability."
- "Understand coding pitfalls that reduce portability."
- "Understand coding pitfalls that reduce maintainability."
keypoints:
- "Nextflow is not sensitive to whitespace. Use it to layout code for readability."
- "Use comments and whitespace to group chunks of code to describe big picture functionality."
- "Avoid `params.parameter` in a process. Pass all parameters using input channels."
- "Input files should be passed using input channels."
---

## Nextflow coding practices

Nextflow is a powerful flexible language that one can code in a variety of ways.
This can lead to poor practices in coding. For example, this can lead
to the workflow only working under certain configurations or execution platforms.
Alternatively, it can make it harder for someone to contribute to a codebase.
These are some useful coding tips that make maintaining and porting your
workflow easier. 


{% include links.md %}
