rule '.md' => '.Rmd' do |t|
  sh "Rscript -e \"library(knitr); knit('#{t.source}', '#{t.name}')\""
end

rule '.pdf' => '.md' do |t|
  sh "pandoc -V geometry:margin=1in --highlight-style tango #{t.source} -o #{t.name}"
end

rule '.docx' => '.md' do |t|
  sh "pandoc #{t.source} -o #{t.name}"
end

