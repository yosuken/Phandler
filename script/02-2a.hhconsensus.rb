
name, ftmp, fcfa = ARGV

if File.exist?(ftmp) #and !(File.exist?(fcfa))
  ### ftmp: align.tmp (multi-line fasta format, including small case letters)
  ### fout: align.w_cons.fa (single-line fasta format, all uppercase)

  open(fcfa, "w"){ |fw|
    open(ftmp){ |fr|
      flag = 0
      seq  = ""
      # fr.gets ### discard the first line begins with "#"

      while l = fr.gets
        if l[0] == ">"
          if flag == 0 ### for the first entry
            gid, *info = l[1..-1].chomp.split(" ")

            # [1] rename consensus seq as #{name}_consensus
            # gid  = "#{name}_consensus"

            # [2] rename consensus seq as #{name}
            gid  = "#{name}"

            lab  = info.size > 0 ? "#{gid} #{info*" "}" : gid
            flag = 1
          else
            lab  = l[1..-1]
          end

          fw.puts seq if seq != ""
          seq = ""

          fw.puts ">#{lab}"
        else
          seq += l.strip.upcase
        end
      end

      fw.puts seq if seq != ""
      seq = ""
    }
  }
end
