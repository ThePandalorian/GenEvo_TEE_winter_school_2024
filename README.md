# GenEvo/TransEvo winter school (2024)

Instructors: Hanna Kokko, Thomas Keaney, Hanna ten Brink and Victor Ronget

The purpose of this course is to introduce mathematical modelling to researchers that primarily conduct empirical research in the fields of ecology, evolution and molecular biology. The course should also be valuable to young theoreticians looking to expand their modelling toolkit. 

This github repository is the place you should look for all online resources associated with the school. Materials associated with each workshop can be downloaded from here, by 'cloning' the repository, or downloading it in ZIP form. To do this, look for the green `code` icon on the top right of this page. Click and copy the web URL to clone, or simply download the ZIP. 

## Required software 

The school is divided into four workshops, each led by a different instructor. These will involve coding in real time, using the `R-studio` IDE to run `R`. To participate in these workshops, you'll therefore need to:

1. Have installed a recent version of R: download it [here](https://cran.rstudio.com/).
2. Have R-studio installed: download it [here](https://posit.co/download/rstudio-desktop/).

`R` is open source, which, among other things, means it's completely free to download. If you have any trouble with this, let Tom know on the Sunday evening and he'll give you a hand. 

## How we'll code 

Just by downloading the files in this repository, you'll be ready to start coding in a transparent, reproducible way, which is fast becoming an essential skill in our fields. We'll work using an `R.project`, which removes any need to muck around with working directories or paths, which is nice because your code becomes user-agnostic. We'll encourage you to write code in `quarto` documents, which are a modern version of R-markdown docs, and can be used simultaneously as an `R.script` and text editor. That means you'll leave each workshop with a polished report that includes all your notes, documents your code, and displays the figures it produces. 

If you're already familiar with workflows like this, that's great! But if this sounds completely foreign, don't worry! That's why we're here. 

## Workshop 1: Discrete Elephants

**Instructor**: Tom Keaney

I will introduce you to modelling by considering an intriguing observation from Gorongosa national park in Mozambique. There, the African Elephant population has shown rapid evolution of tusklessness from 1970 to the mid-2010s. Together, we will build a simple population genetic model to describe the evolution of tusklessness, while making many simplifying assumptions. As we progress through the exercise, we will introduce general modelling guidelines, explain why models are valuable in the first place, and hopefully make clear the strengths and weaknesses of adopting a relatively simple approach. 

## Workshop 2: Continuous Elephants

**Instructor**: Hanna Kokko

On day two, we will take our elephant model and make it more realistic. I will ask you to think about the simplifying assumptions made in Workshop 1, and whether or not they adequately capture the life history of a long-lived organism. We will then explicitly incorporate some life history into the model, and explore how this changes the predictions the model produces.

## Workshop 3: Responding to change

**Instructor**: Hanna ten Brink

All organisms must cope with environmental change, which can be either predictable (e.g., the daily cycle of sunrise and sunset) or unpredictable (e.g., next week's weather). These changes occur on different timescales, from rapid changes (e.g., a rain event) to gradual trends (e.g., seasonal change). First, I will introduce you to the different strategies or response modes that organisms adopt to deal with such changes. We’ll together discuss under which conditions each strategy is adaptive, and when it might not be. Next, we will brainstorm the key components needed to build a model to test these ideas. You’ll then work with an existing model to investigate our hypothesis on the conditions under which each response mode is likely to evolve.  

This model will focus on phenotypic evolution, which is one approach to understand how organisms adapt to change. However, this approach simplifies reality by leaving out the regulatory mechanisms that underlie these evolutionary responses. To address this, I will briefly introduce gene regulatory network models as a way to study evolutionary responses to changing environments. While we won’t have time for hands-on work on this topic, this introduction will give you a start for further exploration on your own.

## Workshop 4: Ageing Dolphins

**Instructor**: Victor Ronget

It is well known that individuals in the same population can exhibit marked heterogeneity in terms of vital rates (this is a demography term for things like birth and death rates). One of the most important factors in this demographic heterogeneity is age. We expect, for example, changes in mortality and fertility rates during juvenile growth or for older senescent individuals. We will first learn how to accurately describe these age-specific changes in mortality and fertility using life table analyses from several mammalian demographic datasets. We will discuss the modeling choices and assumptions associated with these methods. We will then use these age-specific trajectories to build population projection models ranging from discrete to continuous time scales, with particular emphasis on using these models to describe the main evolutionary theories of aging.


