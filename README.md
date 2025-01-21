## Overview

`dnaapler_wrapper.py` is a Python script that respects the `[topology=...]` tags of contig headers:

- `[topology=circular]`: contigs are rotated (default Dnaapler behavior)
- `[topology=linear]`: contigs are reverse-complemented
- No topology tag: behavior depends on `--default-topology`; by default it raises an error

It is inspired by this issue: https://github.com/gbouras13/dnaapler/issues/82

## Installation

Simply run `dnaapler_wrapper.py` in the same environment where `dnaapler` is installed. It requires no additional
dependencies.

## Usage

Try `python dnaapler_wrapper.py --help` to see the available options.

## Reasoning

I use NCBI PGAP FASTA headers after assembly and this may look as follows:

```fasta
>contig_1 [topology=circular]
...
>contig_2 [topology=linear]
...
>contig_3
...
```

Dnaapler does not respect the `[topology=...]` tags and will simply rotate all of these contigs.

I only want it to rotate the contigs that are tagged as circular. This is what `dnaapler_wrapper.py` does.

The output for contigs 1 and 2 will be:

```fasta
>contig_1 [topology=circular] [dnaapler=rotated] [dnaapler-gene=repA]
...
>contig_2 [topology=linear] [dnaapler=reverse-complement] [dnaapler-gene=repA]
...
```

For the ambiguous `contig_3`, the output depends on `--default-topology`:

1) `raise`: `dnaapler_wrapper.py` crashes: `ValueError: Unknown topology in contig_3`
2) `linear`: `>contig_3 [dnaapler=reverse-complement] [dnaapler-gene=repA]`
3) `circular`: `>contig_3 [dnaapler=rotated] [dnaapler-gene=repA]`
