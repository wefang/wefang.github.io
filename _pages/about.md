---
layout: about
title: About
permalink: /
subtitle: >
  Postdoctoral Fellow, Department of Biomedical Engineering, Johns Hopkins University<br><span class="address">855 N. Wolfe St., Rangos 570, Baltimore, MD 21205 · <a href="#" class="al-email-protect" data-eu="wfang9" data-ed="jh.edu"><span class="al-email-text">wfang9</span><span class="al-email-sep"> [at] </span><span class="al-email-text">jh</span><span class="al-email-sep"> [dot] </span><span class="al-email-text">edu</span></a></span>

profile:
  align: right
  image: prof_pic.jpg
  image_circular: true

selected_papers: true
social: true

announcements:
  enabled: false

latest_posts:
  enabled: false
---

**I am on the academic job market in 2026–27**, applying for faculty positions in computational biology, biomedical engineering and biostatistics. My [CV is here](/assets/pdf/Weixiang_Fang_CV.pdf).

I am a postdoctoral fellow in the [Kalhor Lab](https://kalhorlab.bme.jhu.edu) at Johns Hopkins, where I work closely with wet lab scientists. I received my **PhD in Biostatistics** from the Johns Hopkins Bloomberg School of Public Health, working with [Dr. Hongkai Ji](https://jilab.org).

I work with genomic lineage recorders, which can now barcode millions of cells and trace the division tree of a whole embryo. Learning biology from that tree is much harder. I developed **[Quantitative Fate Mapping](/publications/)** (*Cell*, 2022) to do it, recovering when progenitors commit and how many cells go to each fate, decisions that cannot be watched otherwise.

I also built **[FUNCODE](/software/)** for the ENCODE Consortium, because regulatory DNA evolves fast, and an active mouse enhancer often aligns to a human sequence with no function at all. FUNCODE scores whether an element's functional signal is conserved between human and mouse, which helps prioritize mouse findings for translation to human.

### Research Areas

<div class="research-list">
  <div class="research-item">
    <h4>Reading the history of cells from their DNA</h4>
    <p>Engineered recorders now write each cell's history into its own DNA as it divides. I build the models that read those records back into the tree of divisions that built a tissue, at the scale of whole embryos.</p>
    <img src="/assets/img/research/lineage.png" alt="A lineage tree with colored recording events on its branches, and the character matrix of recorded sites in each sampled cell">
  </div>
  <div class="research-item">
    <h4>What decides a cell's fate</h4>
    <p>A cell's future depends both on its own state and on the signals around it. I build models of development that ask how much of a cell's fate is driven by what signals, and how much can still be changed.</p>
    <img src="/assets/img/research/fate.png" alt="A simulation of clones growing across a tissue with a signal gradient, colored by clone and by fate">
  </div>
  <div class="research-item">
    <h4>From model organisms to human</h4>
    <p>Much of what we know about how fate is controlled was learned in animals. I develop methods that carry that knowledge into human, and show where it does not carry.</p>
    <img src="/assets/img/research/transfer.png" alt="A phylogeny of zebrafish, mouse, macaque and human with functional genomics tracks, syntenic regions, and each organism's proposal onto the human track">
  </div>
</div>
