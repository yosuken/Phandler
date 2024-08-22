
fgap, fhhr, fpos = ARGV

# fgap: consensus sequence (gapped)
# fhhr: hhalign result
gapped_seq = IO.readlines(fgap)[1].strip
nogap_seq  = gapped_seq.gsub("-", "")

n2g = {} ### 1-based position convert from nogap to gapped
j = 0 ### nogap_pos
(0..gapped_seq.size-1).each{ |i| ## gapped_pos
  if gapped_seq[i] == "-"
  else
    n2g[j+1] = i+1 ### store 1-based conversion
    j += 1
  end
}

# {{{ [0] range function
class Range
  include Comparable

  def <=>(other)
    self.min <=> other.min
  end
  def overlap?(t)
    s = self
    if (s.max - t.min) * (t.max - s.min) < 0 ## no overlap
      return false
    else
      return true
    end
  end
  def &(t)
    s = self
    if (s.max - t.min) * (t.max - s.min) < 0 ## no overlap
      return nil
    else
      return ([s.min, t.min].max)..([s.max, t.max].min)
    end
  end
  def |(t)
    s = self
    if (s.max - t.min) * (t.max - s.min) < 0 ## no overlap
      return [s, t]
    else
      return [([s.min, t.min].min)..([s.max, t.max].max)]
    end
  end
  def include?(t)
    s = self
    o = s & t  ## overlap of s and t
    return false unless o
    if o.min == t.min and o.max == t.max ## s includes t
      return true
    else
      return false
    end
  end
end
def merge_ranges(ranges) ## ranges = [3..300, 320..500, 504..732, ...]
  return ranges if ranges.size < 2

  rs  = ranges.sort
  (-rs.size..-2).each{ |j| ## index from right
    merged = rs[j] | rs[j+1]
    case merged.size
    when 1 ## overlap detected
      rs[j] = merged[0]
      rs.delete_at(j+1)
    when 2 ## overlap not detected
      ## do nothing
    else raise
    end
  }
  return rs
end
# }}} range function


hhrH = {} ### nogap_pos --> hhr_pos (1-baed)
hhr  = "" ### store hhr_aa
# {{{ parse_hhr(fhhr, hhrH, hhr)
def parse_hhr(fhhr, hhrH, hhr)

  ### query is nogap consensus seq. target is hhr.
  que, tar = "", ""
  que_off, tar_off = 0, 0 ### offset (the number of amino acid befor the alignment)

  flag = 0
  IO.readlines(fhhr).each{ |l|
    if l[0..3] == "No 1"
      flag = 1
    elsif flag == 1 and l[0, 7] == "Probab="
      flag = 2
    elsif flag >= 2 and l[0..2] == "No "
      break ### 2nd entry --> break
    elsif flag >= 2 and l.strip == ""
      flag = 2
    elsif flag >= 2
      flag += 1

      alpha, name, start, seq, stop, len = l.strip.split(/\s+/)
      case flag
      when 3 ### query line
        ### offset (the number of amino acid befor the alignment)
        if que.size == 0
          que_off = start.to_i - 1
        end
        que += seq
      when 7 ### target line
        ### offset (the number of amino acid befor the alignment)
        if tar.size == 0
          tar_off = start.to_i - 1
        end
        tar += seq
      end 
    end
  }

  ### fill hhr
  tar_off.times do hhr << "-" end

  ### make hhrH and hhr
  (0..tar.size-1).each{ |i|
    if tar[i] == "-"
      tar_off -= 1
    end
    if que[i] == "-"
      que_off -= 1
    end

    i2 = i + 1 + tar_off ## 1-based position of target
    j2 = i + 1 + que_off ## 1-based position of query

    if tar[i] != "-" and que[i] != "-"
      hhrH[j2] = i2
    end
    if tar[i] != "-"
      # hhr << que[i]
      hhr << tar[i]
    end
  }

  return 
end
# }}}
parse_hhr(fhhr, hhrH, hhr)
p hhrH

out = []
(1..nogap_seq.size).each{ |i| ## 1-based
  nogap_pos   = i
  gapped_pos  = n2g[i] || "-"
  hhr_pos     = hhrH[i]

  nogap_aa    = nogap_seq[i - 1]
  gapped_aa   = gapped_seq[gapped_pos - 1]
  hhr_aa      = hhr_pos ? hhr[hhr_pos - 1] : "."

  hhr_pos ||= "-" ### not assigned

  out << [gapped_pos, gapped_aa, nogap_pos, nogap_aa, hhr_pos, hhr_aa]*"\t"
}

open(fpos, "w"){ |fw|
  header = %w|aligned_pos aligned_aa nogap_pos nogap_aa hhr_pos hhr_aa|*"\t"
  fw.puts header
  fw.puts out
}
