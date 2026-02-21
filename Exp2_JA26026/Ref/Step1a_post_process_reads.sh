zcat Merged/$1.un.1.fq.gz Merged/$1.fq.gz |awk '{ID=$0;getline;Seq=$0;getline;getline;if (length(Seq)>47 && length(Seq)<53) {print ID"\n"Seq"\n+\n"$0}}' | gzip > Merged/$1.processed.fq.gz
