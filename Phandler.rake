

# {{{ procedures
WriteBatch  = lambda do |t, jobdir, outs|
	jdir = "#{jobdir}/#{t.name.split(":")[-1]}"; mkdir_p jdir unless File.directory?(jdir)
  jnum = outs.size

  if jnum > 0
    outs.each_slice(jnum).with_index(1){ |ls, idx| ## always 1 file
      open("#{jdir}/#{t.name.split(".")[1..-1]*"."}.sh", "w"){ |fjob|
        fjob.puts ls
      }
    }
  else
    open("#{jdir}/#{t.name.split(".")[1..-1]*"."}.sh", "w"){ |fjob| } ## clear job file
  end
end

RunBatch    = lambda do |t, jobdir, ncpu, logdir|
	jdir = "#{jobdir}/#{t.name.split(":")[-1]}"
  ldir = "#{logdir}/#{t.name.split(":")[-1]}"; mkdir_p ldir unless File.directory?(ldir)

  Dir["#{jdir}/*.sh"].sort_by{ |fin| fin.split(".")[-1].to_i }.each{ |fin| ## always 1 or 0 file
    next if File.zero?(fin)
    sh "parallel --jobs #{ncpu} --joblog #{ldir}/parallel.log <#{fin}"
	}
  open("#{ldir}/exit", "w"){ |fw| fw.puts "exit at #{Time.now.strftime("%Y-%m-%d_%H:%M:%S")}" }
end

PrintStatus = lambda do |current, total, status, t|
	puts ""
	puts "\e[1;32m===== #{Time.now}\e[0m"
	puts "\e[1;32m===== step #{current} / #{total} (#{t.name}) -- #{status}\e[0m"
	puts ""
	$stdout.flush
end

CheckVersion = lambda do |commands, cmd2path|
	commands.each{ |command|

    path = $cmd2path[command]||command

    f = if path.split("/").size == 1 then `which #{path}`.strip
        elsif File.exist?(path) then path
        else `which #{path}`.strip
        end

		str = case command
					when "ruby"
						%|#{path} --version 2>&1|
					when "parallel"
						# %{which parallel && LANG=C parallel --version 2>&1 |head -n 1}
						%{LANG=C #{path} --version 2>&1 |head -n 1}
					when "hmmbuild", "hmmsearch"
						%{#{path} -h 2>&1 |head -n 2}
					when "hhalign"
						%{#{path} -h 2>&1 |grep "^HHalign"}
					when "hhsearch"
						%{#{path} 2>&1 |grep "^HHsearch"}
					when "hhblits"
						%{#{path} 2>&1 |grep "^HHblits"}
					when "hhconsensus"
						%{#{path} 2>&1 |grep "^HHconsensus"}
					when "hhmake"
						%{#{path} 2>&1 |grep "^HHmake"}
          when "ffindex_build", "ffindex_apply", "cstranslate", "ffindex_apply_mpi", "cstranslate_mpi", "ffindex_order", "reformat.pl"
						%{#{path} 2>&1 |grep -i "^Usage"}
					when "mmseqs"
						%{#{path} -h |grep "^MMseqs2 Version:"}
					when "weblogo"
						%{#{path} --version}
          when "CPAP"
						%{#{path} -v 2>&1}
					end
		puts ""
		puts "\e[1;32m===== check version: #{command}\e[0m"
		puts ""
		puts "path: #{f}"
		puts ""
		puts "$ #{str}"
		### run
		out = `#{str}`
    if $?.exitstatus == 0 ### zero exit
      puts out
    else ### nonzero exit
      raise("#{Errmsg} #{command}: nonzero exit status #{$?.exitstatus}. command may not be found.")
    end
		### flush
		$stdout.flush
	}
end
# }}} procedures


# {{{ task controller
task :default do
	### constants from arguments
	Odir      = ENV["outdir"]           ## output directory
	OdirExist = ENV["outdir_exist"]     ## output directory exist? ("true" or "")
	Fins      = ENV["input"]            ## input MSA files
  Evalue    = ENV["evalue"].to_f      ## num of CPUs
  Ncpu      = ENV["ncpus"].to_i       ## num of CPUs
  Iter      = ENV["iter"].to_i        ## hhblits iteration (default: 2)

	### constants
  Z           = 1_000_000 # hmmsearch data size for evalue calculation
  Errmsg      = "\e[1;31mError:\e[0m"
  Wrnmsg      = "\e[1;35mWarning:\e[0m"
  MmCov       = "0.5"
  MmMemG      = 20
  MmNcpu      = 1 ### mmseqs
  HmNcpu      = 1 ### hmmsearch
  HhNcpu      = 1 ### hhblits
  MmIdts      = %w|0.9 0.8 0.7 0.6 0.5 0.4 0.3|
  NumStep     = ""
  NseqThre    = 1500 ### threshold for CPAP: if the seq size is <= NseqThre, CPAP will be executed.
  # if NseqThre == 2000, then raised error. ==> "Error: C stack usage  7972496 is too close to the limit"
  ClustMethod = "ALL"
  HhNHit      = 50

	### check version
  cmds      = %w|ruby parallel hmmbuild hmmsearch hhalign hhsearch hhblits hhconsensus hhmake ffindex_build ffindex_apply cstranslate ffindex_order reformat.pl mmseqs CPAP weblogo|
  # cmds      = %w|ruby parallel hmmbuild hmmsearch hhalign hhsearch hhblits hhconsensus hhmake ffindex_build ffindex_apply_mpi cstranslate_mpi ffindex_order reformat.pl CPAP|
  $cmd2path = {}
  cmds.each{ |cmd|
    if cmd == "CPAP"
      $cmd2path[cmd] = "/aptmp/yosuke/dev/tool/CPAP/CPAP"
    else
      $cmd2path[cmd] = cmd 
    end
  }
	CheckVersion.call(cmds, {  })

  ## dir path
  Scrdir    = "#{File.dirname(__FILE__)}/script"
  Jobdir    = "#{Odir}/batch"
  Tmpdir    = "#{Odir}/tmp"
  Logdir    = "#{Odir}/log/tasks"
  Resdir    = "#{Odir}/result"     

  Seqdir    = "#{Odir}/result/seq"     
  SeqAdir   = "#{Odir}/result/seq_all"     
  Clstdir   = "#{Odir}/result/cluster"     
  ClstAdir  = "#{Odir}/result/cluster_all"
  Heatdir   = "#{Odir}/result/heatmap"     
  HeatAdir  = "#{Odir}/result/heatmap_all"     
  Logodir   = "#{Odir}/result/seqlogo"     

  Hmmdir    = "#{Odir}/result/hmm"
  HmmAdir   = "#{Odir}/result/hmm_all"
  HmmA      = "#{Odir}/result/hmm_all/all.hmm"
  HmmSdir   = "#{Odir}/result/hmmsearch" ### hmmsearch result

  Hhdbdir   = "#{Odir}/result/hhdb"
  HhdbAdir  = "#{Odir}/result/hhdb_all"
  Hhdb      = "#{Odir}/result/hhdb_all/all"
  Hhmdir    = "#{Odir}/result/hhm"
  Hhbdir    = "#{Odir}/result/hhblits"
  Hhsdir    = "#{Odir}/result/hhsearch"
  Hhadir    = "#{Odir}/result/hhalign"

  ## validate evalue
  raise("#{Errmsg} -e, --evalue '#{ENV["evalue"]}' should be positive number.") if Evalue <= 0

  ## Odir exist?
  $stderr.puts "\n\n#{Wrnmsg} output directory #{Odir} already exists. Overwrite it.\n\n" if OdirExist == "true"

	### define tasks
  $tasks, $fpara = {}, {} ### tasks, parallel_flag
  $tasks["TOP"] , $fpara["TOP"]  = %w|T1 T2|, false
  $tasks["T1"]  , $fpara["T1"]   = %w|01-1a.validate_input|, false
  $tasks["T2"]  , $fpara["T2"]   = %w|T2-1 T2-2 T2-3 T2-4|, true
  $tasks["T2-1"], $fpara["T2-1"] = %w|02-1a.hmmbuild 02-1b.concat_hmm 02-1c.hmmsearch|, false ### hmmer
  $tasks["T2-2"], $fpara["T2-2"] = %w|02-2a.hhconsensus 02-2b.hhmake 02-2c.hhdb 02-2d.concat_hhdb 02-2e.hhblits 02-2f.hhsearch 02-2g.parse_hhblits_hhsearch 02-2h.cons_hhalign|, false ### hh-suite

  $tasks["T2-3"], $fpara["T2-3"] = %w|02-3a.mmseqs_cluster 02-3b.CPAP|, false ### heatmap
  # $tasks["T2-3"], $fpara["T2-3"] = %w|02-3a.mmseqs_cluster|, false ### heatmap

  $tasks["T2-4"], $fpara["T2-4"] = %w|02-4a.weblogo 02-4b.skylign|, false ### seqlogo

  def do_task(t, n)
    if Rake::Task.task_defined?(t)
      Rake::Task[t].invoke(t, n)
    else
      if $tasks[t]
        if $fpara[t]  ### parallel
          # Parallel.each($tasks[t], in_processes: [$tasks[t].size, Ncpu].min){ |_t|
          #   do_task(_t)
          # }
          $tasks[t].each.with_index(1){ |_t, idx|
            do_task(_t, n)
          }
        else
          $tasks[t].each.with_index(1){ |_t, idx|
            do_task(_t, n)
          }
        end
      end
    end
  end
   
  ### run
  do_task("TOP", Ncpu)

end
# }}} default (run all tasks)


# {{{ function
# {{{ cmd
def cmd(cmd)
  path = $cmd2path[cmd]

  raise("#{cmd}: path of the command is not defined.") unless path

  return path
end
# }}}
# }}}


# {{{ tasks
# {{{ 01-1a.validate_input 
desc "01-1a.validate_input"
task "01-1a.validate_input", ["step", "ncpu"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)

  ### expand tilde and wildcard
  fins = Fins.split(",").inject([]){ |a, path| a += Dir[path.gsub("~", ENV["HOME"])].sort }

  fins = fins.map{ |fin|
    if File.zero?(fin)
      $stderr.puts "#{Wrnmsg} #{fin} is empty. Skip this file."
      next nil
    else
      File.open(fin){ |fr|
        if fr.gets[0] != ">"
          $stderr.puts "#{Wrnmsg} #{fin} is not valid fasta format. The first character should be '>'. Skip this file."
          next nil
        end
      } 
    end
    fin ## valid fasta
  }.compact

  $stderr.puts ["", "", "\e[1;32m===== query file validated (valid file: n=#{fins.size}) \e[0m"]
  raise("#{Errmsg} no query file detected.") if fins.size == 0

  ### store file info
  $input = []
  names  = {}

  fins.each.with_index(1){ |fin, idx|
    ### file names should have valid extension.
    fname = File.basename(fin) 
    if fname.split(".").size <= 1
      raise("'#{fname}': file name should have any extension such as .fa, .fasta, and .mfa.") if fins.size == 0
    end

    ### file names should be uniq.
    name  = fname.split(".")[0..-2]*"."
    if names[name]
      raise("'#{name}' is used as a name of the given path #{fin}, but the name found twice within input files. The name should be uniq.")
    end
    names[name] = 1

    ### open file and parse seq length
    gid2len = {}
    open(fin){ |fr|
      gid, len = "", 0
      while l = fr.gets 
        if l[0] == ">"
          gid2len[gid] = len if gid != ""

          gid = l.strip.split(" ")[0][1..-1]
          len = 0
        else
          len += l.strip.size
        end
      end
      gid2len[gid] = len if gid != ""
      len = 0
    }
    # require 'pp' ; pp gid2len

    ### input is alignment, so length should be uniform
    if lens = gid2len.values.uniq and lens.size > 1
      raise("'#{fname}': the file includes sequences of different lengths. An MSA file should have sequences of the same length.")
    end
    len   = lens[0]
    n_seq = gid2len.size

    ### open file and write no gapped fasta
    sdir  = "#{Seqdir}/#{name}" ; mkdir_p sdir unless File.directory?(sdir)
    faln  = "#{sdir}/align.fa"
    fa    = "#{sdir}/nogap.fa"
    fa3m  = "#{sdir}/align.a3m"
    fcfa  = "#{sdir}/align.w_cons.fa" ### alinged fasta with consensus seq at the top.
    fhmm  = "#{Hmmdir}/#{name}.hmm"
    fhhm  = "#{Hhmdir}/#{name}.hhm"
    fhhrS = "#{Hhsdir}/raw/#{name}.hhr"
    fhhrB = "#{Hhbdir}/raw/#{name}.hhr"

    fw0 = open(faln, "w")
    fw1 = open(fa  , "w")
    open(fin){ |fr|
      seq = ""
      while l = fr.gets 
        if l[0] == ">"
          fw0.puts seq if seq != ""
          fw1.puts seq.gsub("-", "") if seq != ""
          seq = ""

          fw0.puts l
          fw1.puts l
        else
          seq += l.strip
        end
      end
      fw0.puts seq if seq != ""
      fw1.puts seq.gsub("-", "") if seq != ""
      seq = ""
    }
    [fw0, fw1].each{ |fw| fw.close }


    h = { idx: idx, name: name, faln: faln, fa: fa, fa3m: fa3m, fcfa: fcfa, fhmm: fhmm, fhhm: fhhm, fhhrS: fhhrS, fhhrB: fhhrB, length: len, n_seq: n_seq }
    $stderr.puts h.inspect
    $input << h
  }
end
# }}}

# {{{ 02-1a.hmmbuild
desc "02-1a.hmmbuild"
task "02-1a.hmmbuild", ["step", "ncpu"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	outs   = []

  mkdir_p Hmmdir unless File.directory?(Hmmdir)

  $input.each{ |set|
    name = set[:name]
    faln = set[:faln]
    fhmm = set[:fhmm]
    flog = "#{fhmm}.log"

    outs << "#{cmd("hmmbuild")} -n #{name} #{fhmm} #{faln} >#{flog} 2>&1" unless File.exist?(fhmm)
  }

	WriteBatch.call(t, Jobdir, outs)
	RunBatch.call(t, Jobdir, args.ncpu, Logdir)
end
# }}}

# {{{ 02-1b.concat_hmm
desc "02-1b.concat_hmm"
task "02-1b.concat_hmm", ["step", "ncpu"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)

  mkdir_p HmmAdir unless File.directory?(HmmAdir)
  next if File.exist?(HmmA)

  open(HmmA, "w"){ |fw|
    $input.each{ |set|
      fhmm = set[:fhmm]

      fw.puts IO.read(fhmm)
    }
  }
end
# }}}

# {{{ 02-1c.hmmsearch
desc "02-1c.hmmsearch"
task "02-1c.hmmsearch", ["step", "ncpu"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	outs   = []

  $input.map{ |set|
    name = set[:name]
    fhmm = set[:fhmm]
    fa   = set[:fa]

    odir = "#{HmmSdir}/#{name}" ; mkdir_p odir unless File.directory?(odir)
    fout = "#{odir}/hmmsearch.out"
    flog = "#{odir}/hmmsearch.log"

    fgz  = "#{fout}.gz"
    next if File.exist?(fgz)

    outs << "#{cmd("hmmsearch")} --cpu #{HmNcpu} -Z #{Z} -o #{fout} #{fhmm} #{fa} >#{flog} 2>&1 && gzip #{fout}"
  }

	WriteBatch.call(t, Jobdir, outs)
	RunBatch.call(t, Jobdir, args.ncpu, Logdir)
end
# }}}

# {{{ 02-2a.hhconsensus
desc "02-2a.hhconsensus"
task "02-2a.hhconsensus", ["step", "ncpu"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	outs   = []
  script = "#{Scrdir}/#{t.name}.rb"

  $input.each{ |set|
    name  = set[:name]
    faln  = set[:faln]
    fa3m  = set[:fa3m]
    fcfa  = set[:fcfa]
    ftmp  = fa3m.gsub(/\.a3m$/, ".tmp")

    next if File.exist?(fcfa)

    out   = []
    out  << "#{cmd("hhconsensus")} -M 50 -maxseq 65535 -maxres 65535 -i #{faln} -o #{fa3m} >#{fa3m}.log 2>&1"
    out  << "#{cmd("reformat.pl")} a3m fas #{fa3m} #{ftmp} >#{ftmp}.log 2>&1"
    out  << "#{cmd("ruby")} #{script} #{name} #{ftmp} #{fcfa} >#{fcfa}.log 2>&1" ### rename the consensus entry (the 1st entry)

    outs << out*" && "
  }

	WriteBatch.call(t, Jobdir, outs)
	RunBatch.call(t, Jobdir, args.ncpu, Logdir)
end
# }}}

# {{{ 02-2b.hhmake
desc "02-2b.hhmake"
task "02-2b.hhmake", ["step", "ncpu"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	outs   = []

  mkdir_p Hhmdir unless File.directory?(Hhmdir)

  $input.each{ |set|
    name  = set[:name]
    fhhm  = set[:fhhm]
    fcfa  = set[:fcfa]

    flog  = "#{fhhm}.log"

    outs << "#{cmd("hhmake")} -maxres 65535 -i #{fcfa} -o #{fhhm} -M 50 >#{flog} 2>&1" unless File.exist?(fhhm)
  }

	WriteBatch.call(t, Jobdir, outs)
	RunBatch.call(t, Jobdir, args.ncpu, Logdir)
end
# }}}

# {{{ 02-2c.hhdb
desc "02-2c.hhdb"
task "02-2c.hhdb", ["step", "ncpu"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	outs   = []

  $input.each{ |set|
    name  = set[:name]
    fhhm  = set[:fhhm]
    fcfa  = set[:fcfa]

    odir  = "#{Hhdbdir}/#{name}"
    mdir  = "#{Hhdbdir}/#{name}/msa"
    next if File.exist?("#{odir}/#{name}_a3m.ffdata")

    out  = []
    out << "mkdir -p #{mdir}"
    out << "cp -p #{fcfa} #{mdir}/#{name}.fa" ### rename: align.w_cons.fa --> #{name}.fa
    out << "( pushd #{mdir} >/dev/null"
    out << "#{cmd("ffindex_build")} -s ../#{name}_msa.ff{data,index} ."
    out << "cd .."
    out << "#{cmd("ffindex_apply")} #{name}_msa.ff{data,index} -i #{name}_a3m_wo_ss.ffindex -d #{name}_a3m_wo_ss.ffdata -- #{cmd("hhconsensus")} -M 50 -maxres 65535 -i stdin -oa3m stdout -v 0"
    out << "mv #{name}_a3m_wo_ss.ffdata #{name}_a3m.ffdata"
    out << "mv #{name}_a3m_wo_ss.ffindex #{name}_a3m.ffindex"
    out << "#{cmd("ffindex_apply")} #{name}_a3m.ff{data,index} -i #{name}_hhm.ffindex -d #{name}_hhm.ffdata -- #{cmd("hhmake")} -i stdin -o stdout -v 0"
    out << "#{cmd("cstranslate")} -f -x 0.3 -c 4 -I a3m -i #{name}_a3m -o #{name}_cs219"

    out << "sort -k3 -n -r #{name}_cs219.ffindex | cut -f1 > sorting.dat"
    out << "#{cmd("ffindex_order")} sorting.dat #{name}_hhm.ff{data,index} #{name}_hhm_ordered.ff{data,index}"
    out << "mv #{name}_hhm_ordered.ffindex #{name}_hhm.ffindex"
    out << "mv #{name}_hhm_ordered.ffdata #{name}_hhm.ffdata"

    out << "#{cmd("ffindex_order")} sorting.dat #{name}_a3m.ff{data,index} #{name}_a3m_ordered.ff{data,index}"
    out << "mv #{name}_a3m_ordered.ffindex #{name}_a3m.ffindex"
    out << "mv #{name}_a3m_ordered.ffdata #{name}_a3m.ffdata"
    out << "popd >/dev/null )"

    outs << out*" && "
  }

	WriteBatch.call(t, Jobdir, outs)
	RunBatch.call(t, Jobdir, args.ncpu, Logdir)
end
# }}}

# {{{ 02-2d.concat_hhdb
desc "02-2d.concat_hhdb"
task "02-2d.concat_hhdb", ["step", "ncpu"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	outs   = []


  odir  = "#{HhdbAdir}"
  mdir  = "#{HhdbAdir}/msa"
  name  = "all"
  next if File.exist?("#{odir}/#{name}_a3m.ffdata")

  out = []
  out << "mkdir -p #{mdir}"
  $input.each{ |set|
    fcfa = set[:fcfa]
    out  << "cp -p #{fcfa} #{mdir}/#{set[:name]}.fa" ### rename: align.w_cons.fa --> #{name}.fa
  }
  out << "( pushd #{mdir} >/dev/null"
  out << "#{cmd("ffindex_build")} -s ../#{name}_msa.ff{data,index} ."
  out << "cd .."
  out << "#{cmd("ffindex_apply")} #{name}_msa.ff{data,index} -i #{name}_a3m_wo_ss.ffindex -d #{name}_a3m_wo_ss.ffdata -- #{cmd("hhconsensus")} -M 50 -maxres 65535 -i stdin -oa3m stdout -v 0"
  out << "mv #{name}_a3m_wo_ss.ffdata #{name}_a3m.ffdata"
  out << "mv #{name}_a3m_wo_ss.ffindex #{name}_a3m.ffindex"
  out << "#{cmd("ffindex_apply")} #{name}_a3m.ff{data,index} -i #{name}_hhm.ffindex -d #{name}_hhm.ffdata -- #{cmd("hhmake")} -i stdin -o stdout -v 0"
  out << "#{cmd("cstranslate")} -f -x 0.3 -c 4 -I a3m -i #{name}_a3m -o #{name}_cs219"

  out << "sort -k3 -n -r #{name}_cs219.ffindex | cut -f1 > sorting.dat"
  out << "#{cmd("ffindex_order")} sorting.dat #{name}_hhm.ff{data,index} #{name}_hhm_ordered.ff{data,index}"
  out << "mv #{name}_hhm_ordered.ffindex #{name}_hhm.ffindex"
  out << "mv #{name}_hhm_ordered.ffdata #{name}_hhm.ffdata"

  out << "#{cmd("ffindex_order")} sorting.dat #{name}_a3m.ff{data,index} #{name}_a3m_ordered.ff{data,index}"
  out << "mv #{name}_a3m_ordered.ffindex #{name}_a3m.ffindex"
  out << "mv #{name}_a3m_ordered.ffdata #{name}_a3m.ffdata"
  out << "popd >/dev/null )"

  outs << out*" && "

	WriteBatch.call(t, Jobdir, outs)
	RunBatch.call(t, Jobdir, args.ncpu, Logdir)
end
# }}}

# {{{ 02-2e.hhblits
desc "02-2e.hhblits"
task "02-2e.hhblits", ["step", "ncpu"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	outs   = []

  opt    = "-mact 0.2" ### default: -mact 0.35

  $input.each{ |set|
    fhhm  = set[:fhhm]
    fhhrB = set[:fhhrB]

    odir = "#{Hhbdir}/raw" ; mkdir_p odir unless File.directory?(odir)

    fout = fhhrB
    flog = "#{fhhrB}.log"
    ferr = "#{fhhrB}.err"
    next if File.exist?(fout)

    outs << "#{cmd("hhblits")} -cpu #{HhNcpu} -n #{Iter} #{opt} -z #{HhNHit} -b #{HhNHit} -i #{fhhm} -d #{Hhdb} -o #{fout} >#{flog} 2>#{ferr}"
  }

	WriteBatch.call(t, Jobdir, outs)
	RunBatch.call(t, Jobdir, args.ncpu, Logdir)
end
# }}}

# {{{ 02-2f.hhsearch
desc "02-2f.hhsearch"
task "02-2f.hhsearch", ["step", "ncpu"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	outs   = []

  opt    = "-mact 0.2" ### default: -mact 0.35

  $input.each{ |set|
    fhhm  = set[:fhhm]
    fhhrS = set[:fhhrS]

    odir = "#{Hhsdir}/raw" ; mkdir_p odir unless File.directory?(odir)

    fout = fhhrS
    flog = "#{fhhrS}.log"
    ferr = "#{fhhrS}.err"
    next if File.exist?(fout)

    outs << "#{cmd("hhsearch")} -cpu #{HhNcpu} #{opt} -z #{HhNHit} -b #{HhNHit} -i #{fhhm} -d #{Hhdb} -o #{fout} >#{flog} 2>#{ferr}"
  }

	WriteBatch.call(t, Jobdir, outs)
	RunBatch.call(t, Jobdir, args.ncpu, Logdir)
end
# }}}

# {{{ 02-2g.parse_hhblits_hhsearch
desc "02-2g.parse_hhblits_hhsearch"
task "02-2g.parse_hhblits_hhsearch", ["step", "ncpu"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)

  def parse_hhr(fhhr, col)
    #  No Hit                             Prob E-value P-value  Score    SS Cols Query HMM  Template HMM
    #   1 IncP-9                         100.0 2.2E-51   2E-52  305.0   0.0  182    2-183     2-183 (183)
    #   2 IncP-10                         90.7 1.3E-05 1.2E-06   45.5   0.0   69   31-99     80-168 (271)
    #   3 IncP-1                          81.8 0.00012 1.1E-05   40.0   0.0   55   25-82    210-271 (302)
    #  ...
    #  20 IncP-10                          0.2       4    0.37   11.9   0.0   14  169-182   109-122 (271)
    #
    # No 1
    # >IncP-9
    h = {}
    flag = 0
    IO.readlines(fhhr).each{ |l|
      if l =~ /^ No Hit /
        flag = 1
      elsif flag == 1 and l.strip =~ /^\d+ /
        name, val = l.strip.split(/\s+/).values_at(1, col)

        h[name] = val unless h[name] ### takes only best hit
      end

      if l[0, 4] == "No 1"
        ### break the loop since each hit report begins
        break
      end
    }

    return h
  end

  names  = $input.map{ |set| set[:name] }

  [:fhhrB, :fhhrS].zip([Hhbdir, Hhsdir]){ |fhhr, odir|
    prefs = %w|probability query_hmm target_hmm evalue score|
    cols  = [2, 8, 9, 3, 5]

    prefs.zip(cols){ |pref, col|
      fout = "#{odir}/#{pref}.tsv"

      data  = {}
      $input.each{ |set|
        next unless File.exist?(set[fhhr])

        h = parse_hhr(set[fhhr], col) ### { <name> => <prob> }
        data[set[:name]] = h
      }

      open(fout, "w"){ |fw|
        fw.puts ["", names]*"\t"
        names.each{ |name|
          o = [name]
          h = data[name]
          names.each{ |_name|
            o << h[_name]||""
          }
          fw.puts o*"\t"
        }
      }
    }
  }
end
# }}}

# {{{ 02-2h.cons_hhalign
desc "02-2h.cons_hhalign"
task "02-2h.cons_hhalign", ["step", "ncpu"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	outs   = []
  script = "#{Scrdir}/#{t.name}.rb"

  names  = $input.map{ |set| set[:name] }
  $input.each{ |set|
    name = set[:name]
    fcfa = set[:fcfa]
    fhhm = set[:fhhm]

    odir = "#{Hhadir}/#{name}"; mkdir_p odir unless File.directory?(odir)

    ### write consensus seq file
    fgap = "#{odir}/cons_gapped.fa"
    fnog = "#{odir}/cons_nogap.fa"
    open(fgap, "w"){ |fw| fw.puts IO.readlines(fcfa)[0..1] }
    open(fnog, "w"){ |fw| fw.puts IO.readlines(fcfa)[0] ; fw.puts IO.readlines(fcfa)[1].gsub("-", "") }

    ### hhalign to self hhm
    fhhr = "#{odir}/cons_nogap.hhr"
    flog = "#{fhhr}.log"
    ferr = "#{fhhr}.err"
    fpos = "#{odir}/cons_position_conv.tsv"

    out = []
    out << "#{cmd("hhalign")} -i #{fnog} -t #{fhhm} -o #{fhhr} >#{flog} 2>#{ferr}"
    out << "#{cmd("ruby")} #{script} #{fgap} #{fhhr} #{fpos}"

    unless File.exist?(fpos)
      outs << out*" && "
    end
  }

	WriteBatch.call(t, Jobdir, outs)
	RunBatch.call(t, Jobdir, args.ncpu, Logdir)
end
# }}}

# {{{ 02-3a.mmseqs_cluster
desc "02-3a.mmseqs_cluster"
task "02-3a.mmseqs_cluster", ["step", "ncpu"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
  outs   = []
  opts   = "--cluster-mode 2 -c #{MmCov} --cov-mode 0 --split-memory-limit #{MmMemG}G --threads #{MmNcpu}" 

  ###
  ### for all
  ###
  name = "all"
  sdir = "#{SeqAdir}/#{name}"; mkdir_p sdir unless File.directory?(sdir)
  fa   = "#{SeqAdir}/#{name}/nogap.fa"

  def add_prefix_and_write(name, fr, fw)
    while l = fr.gets
      if l[0] == ">"
        ### add #{name}_ as prefix of gene id
        fw.puts ">#{name}_#{l[1..-1]}"
      else
        fw.puts l
      end
    end
  end

  if !(File.exist?(fa))
    open(fa, "w"){ |fw|
      $input.each{ |set|
        _fa   = set[:fa]
        _name = set[:name]

        open(_fa){ |fr|
          add_prefix_and_write(_name, fr, fw)
        }
      }
    }
  end

  MmIdts.each{ |idt|
    odir = "#{ClstAdir}/#{name}/idt#{idt}_cov#{MmCov}"; mkdir_p odir unless File.directory?(odir)
    fclu = "#{odir}/out_cluster.tsv"
    next if File.exist?(fclu)

    outs << "pushd #{odir} >/dev/null && ln -s #{File.absolute_path(fa)} && #{cmd("mmseqs")} easy-cluster #{File.basename(fa)} out tmp --min-seq-id #{idt} #{opts} >run.log 2>&1 ; popd >/dev/null"
  }

  ###
  ### for each
  ###
  $input.each{ |set|
    fa   = set[:fa]
    name = set[:name]

    MmIdts.each{ |idt|
      odir = "#{Clstdir}/#{name}/idt#{idt}_cov#{MmCov}"; mkdir_p odir unless File.directory?(odir)
      fclu = "#{odir}/out_cluster.tsv"
      next if File.exist?(fclu)

      outs << "pushd #{odir} >/dev/null && ln -s #{File.absolute_path(fa)} && #{cmd("mmseqs")} easy-cluster #{File.basename(fa)} out tmp --min-seq-id #{idt} #{opts} >run.log 2>&1 ; popd >/dev/null"
    }
  }

	WriteBatch.call(t, Jobdir, outs)
	RunBatch.call(t, Jobdir, args.ncpu, Logdir)
end
# }}}

# {{{ 02-3b.CPAP
desc "02-3b.CPAP"
task "02-3b.CPAP", ["step", "ncpu"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
  outs   = []
  header = %w|name cluster n_sequence used_for_heatmap|
  opts   = "--ncpus 1 --clust-method #{ClustMethod}"

  def count_entries(fa)
    count = 0
    open(fa){ |fr|
      while l = fr.gets
        if l[0] == ">"
          count += 1
        end
      end
    }
    return count
  end

  ###
  ### parse cluster size for all
  ###
  fout = "#{ClstAdir}/num_cluster.tsv"
  # unless File.exist?(fout)
  out   = []
  name  = "all"
  fa    = "#{SeqAdir}/#{name}/nogap.fa"
  n_seq = count_entries(fa)
  fin   = nil

  ### without cluster
  flag = 0
  if n_seq <= NseqThre and flag == 0 and n_seq > 1 ### "y" is just once
    fin = fa
    out << [name, "input", n_seq, "y"]
    flag = 1
  else
    out << [name, "input", n_seq, "n"]
  end

  MmIdts.each{ |idt|
    odir = "#{ClstAdir}/#{name}/idt#{idt}_cov#{MmCov}"; mkdir_p odir unless File.directory?(odir)
    fclu = "#{odir}/out_cluster.tsv"
    frpr = "#{odir}/out_rep_seq.fasta"

    n_seq = IO.readlines(fclu).map{ |l| l.chomp.split("\t")[0] }.uniq.size

    if n_seq <= NseqThre and flag == 0 and n_seq > 1 ### "y" is just once
      fin = frpr
      out << [name, idt, n_seq, "y"]
      flag = 1
    else
      out << [name, idt, n_seq, "n"]
    end
  }

  ### select cluster used for CPAP
  if fin and !(File.directory?("#{HeatAdir}/#{name}"))
    mkdir_p HeatAdir unless File.directory?(HeatAdir)
    outs << "#{cmd("CPAP")} #{opts} #{fin} #{HeatAdir}/#{name} >#{HeatAdir}/#{name}.log 2>&1"
  end

  ### write list
  open(fout, "w"){ |fw|
    fw.puts header*"\t"
    out.each{ |o| fw.puts o*"\t" }
  }
  # end

  ###
  ### parse cluster size for each
  ###
  fout = "#{Clstdir}/num_cluster.tsv"
  # unless File.exist?(fout)
  out   = []

  $input.each{ |set|
    name    = set[:name]
    n_seq   = set[:n_seq]
    fin     = nil

    ### without cluster
    flag = 0
    if n_seq <= NseqThre and flag == 0 and n_seq > 1 ### "y" is just once
      fin = set[:fa]
      out << [name, "input", n_seq, "y"]
      flag = 1
    else
      out << [name, "input", n_seq, "n"]
    end

    MmIdts.each{ |idt|
      odir = "#{Clstdir}/#{name}/idt#{idt}_cov#{MmCov}"; mkdir_p odir unless File.directory?(odir)
      fclu = "#{odir}/out_cluster.tsv"
      frpr = "#{odir}/out_rep_seq.fasta"

      n_seq = IO.readlines(fclu).map{ |l| l.chomp.split("\t")[0] }.uniq.size

      if n_seq <= NseqThre and flag == 0 and n_seq > 1 ### "y" is just once
        fin = frpr
        out << [name, idt, n_seq, "y"]
      else
        out << [name, idt, n_seq, "n"]
      end
    }

    ### select cluster used for CPAP
    if fin and !(File.directory?("#{Heatdir}/#{name}"))
      mkdir_p Heatdir unless File.directory?(Heatdir)
      outs << "#{cmd("CPAP")} #{opts} #{fin} #{Heatdir}/#{name} >#{Heatdir}/#{name}.log 2>&1"
    end
  }

  ### write list
  open(fout, "w"){ |fw|
    fw.puts header*"\t"
    out.each{ |o| fw.puts o*"\t" }
  }
  # end

  ### run batch
	WriteBatch.call(t, Jobdir, outs)
	RunBatch.call(t, Jobdir, args.ncpu, Logdir)
end
# }}}

# {{{ 02-4a.weblogo
desc "02-4a.weblogo"
task "02-4a.weblogo", ["step", "ncpu"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	outs   = []

  $input.each{ |set|
    name  = set[:name]
    faln  = set[:faln]

    odir  = "#{Logodir}/weblogo"; mkdir_p odir unless File.directory?(odir)
    feps  = "#{odir}/#{name}.eps"

    outs << "#{cmd("weblogo")} -n 100 -f #{faln} -F eps >#{feps}" unless File.exist?(feps)
  }

	WriteBatch.call(t, Jobdir, outs)
	RunBatch.call(t, Jobdir, args.ncpu, Logdir)
end
# }}}

# {{{ 02-4b.skylign
desc "02-4b.skylign"
task "02-4b.skylign", ["step", "ncpu"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	outs   = []

  $input.each{ |set|
  }

	WriteBatch.call(t, Jobdir, outs)
	RunBatch.call(t, Jobdir, args.ncpu, Logdir)
end
# }}}
# }}}
