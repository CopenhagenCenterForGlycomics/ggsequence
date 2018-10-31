# ggsequence


Perform sequence alignment and render the results as a ggplot

```R
library(ggplot)
library(ggsequence)
```

## Preparing the alignment

We need to make sure the alignment that we want to plot comes in the right format, which
is as a `MsaAAMultipleAlignment`, with the original sequences as an attribute on the
object

```R
do_alignment = function(sequences) {
	# Make sure sequences have a name
	if (is.null(names(sequences))) {
		names(sequences) = as.character(1:length(sequences))
	}
	# Convert each of the sequence strings to a Biostrings::Views
	seqs = lapply( 1:length(sequences), function(seq_id) {
		views = Biostrings::Views(sequences[seq_id],start=1,end=nchar(sequences[seq_id]))
		names(views) = names(sequences)[seq_id]
		views
	})
	names(seqs) = names(sequences)

	# Perform the MSA on the sequences

	aln = msa::msaClustalOmega( Reduce(c,Map(function(view) { as(view,'AAStringSet')}, seqs)) )

	# Assign the original Views object to the sequences attribute
	attributes(aln)$sequences = seqs
	aln
}
```

Now, we can perform an alignment between two sequences, and plot them

```R
p = ggplot(do_alignment(c(foo='MNTTTMMMNPPPPMNTTTMMMNPPPPMNTTTMMMNPPPPMNTTTMMMNPPPP',bar='NNSMMMPPNNSMMMPPNNSMMMPPNNSMMMPP')))
p
```

Which won't really give us anything useful apart from labels for the sequences and amino acid position. We can add in the sequence.

```R
p + geom_text(aes(label=..aa..))
```

Or, we could add in the sequence, but colour the sequence based upon the conservation of the amino acid using the `conservation` stat.

```R
p + geom_text(aes(label=..aa..,color=..conservation..,),stat="conservation")
```

Or, if we want a quick ClustalW style rendering of sequence

```R
p + geom_text(aes(label=..aa..)) + geom_barcode(overlay=F)
```

Which we could actually put on top of the sequence

```R
p + geom_barcode()
```

Lets say we are only interested in the sequence that is aligned with another sequence - we can use the `gappedSequence` stat to draw segments

```R
p + geom_text(aes(label=..aa..)) + geom_segment(aes(x=..seqstart..,xend=..seqend..),stat="gappedSequence",size=2,colour="black",position = position_nudge(y=-0.1))
```

This is fine, but what if we want to map some data onto our aligned sequences? You can map on a data frame to the alignment using the `alignedSite` stat.
You tell the stat where the data is (`annotations` parameter), and which column the sequence ids are in (`id.column` parameter), and which columns you want available as amino acid positions that have been rescaled
to match the alignment - in this case `start` and `end` (`columns` parameter).

```R
signalpeptide_data=data.frame(seq.ids=c('bar','foo'),start=c(1,1),end=c(3,4))
p + geom_text(aes(label=..aa..)) + geom_segment(aes(x=..start..,xend=..end..),stat="alignedSite",colour="red",size=4,alpha=0.5,annotations=signalpeptide_data,id.column='seq.ids',columns=c('start','end'))
```

Or, you can use some fancy brackets to indicate a region
```R
p + geom_text(aes(label=..aa..)) + geom_bracket(aes(x=..start..,xend=..end..),offset=1,size=0.5,stat="alignedSite",annotations=signalpeptide_data,id.column='seq.ids',columns=c('start','end'))
```

If you have a data frame with specific annotations to add on, you can use the `alignedSite` stat to map that on to the sequences too
```R
site_data = data.frame(seq.ids=c('bar','foo'),site=c(4,5))
site_data$site_label = as.character(site_data$site)
p + geom_text(aes(label=..aa..)) + geom_point(aes(x=..site..),stat="alignedSite",annotations=site_data,id.column='seq.ids',columns=c('site'),position=position_nudge(y=0.1),color='red')
```

Extra columns get automatically made available within the geom, so if you have a column with labels for each row, you can easily access them.
```R
p + geom_text(aes(label=..aa..)) + geom_label(aes(x=..site..,label=..site_label..),stat="alignedSite",annotations=site_data,id.column='seq.ids',columns=c('site'),position=position_nudge(y=0.1),color='red')
```

Or, if you are being really fancy
```R
p + geom_text(aes(label=..aa..)) + ggsugar::geom_sugar(aes(x=..site..),stat="alignedSite",annotations=site_data,id.column='seq.ids',columns=c('site'),position=position_nudge(y=0.1),class='gal(b1-3)galnac')
```

Or, if you are being *EVEN* fancier
```R
site_data = data.frame(seq.ids=c('bar','foo'),site=c(4,5),sugarseq=c('man','gal(b1-3)galnac'))
p + geom_text(aes(label=..aa..)) + ggsugar::geom_sugar(aes(x=..site..,class=..sugarseq..),stat="alignedSite",annotations=site_data,id.column='seq.ids',columns=c('site'),position=position_nudge(y=0.1))
```
