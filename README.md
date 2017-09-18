# RayTracer

A ray tracer application created for the computer graphics course at university.

## Introduction

The code was originally written as my computer graphics homework. Since then, I added new features when I had the time. However, at this state, the code is monolithic, chaotic, and most of the identifiers and comments are written in Hungarian.

## The plan

My plan is to improve the code by

* eliminating the Hungarian identifiers and introduce English ones
* organizing the code into separate files based on the an OO class structure
* finishing unfinishied features (for example, the Hearth shape that has been added is not ready yet)
* fixing bugs in existing features
* refactoring the process of ray tracing:
  * intersect calculation can be simplified to calculate only for unit size of a each object type
  * rays can be transformed into the coordinate space of the model based on its size, rotation, etc