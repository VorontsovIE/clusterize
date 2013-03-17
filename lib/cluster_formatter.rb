class ClusterFormatter
  attr_accessor :clusterer, :branch_len_meth
  def initialize(clusterer, branch_len_meth)
    @clusterer = clusterer
    @branch_len_meth = branch_len_meth
  end
  def compact_name(name)
    mtch = /(?<rep>.+).+?\+\k<rep>/.match(name)
    mtch ? name.gsub(/\+#{mtch[:rep]}/,'+') : name
  end
  
end