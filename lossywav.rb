class Lossywav < Formula
  desc "POSIX port of lossyWAV - a near lossless audio processor."
  homepage "https://github.com/MoSal/lossywav-for-posix"
  head "https://github.com/MoSal/lossywav-for-posix.git"

  depends_on "pkg-config" => :build
  depends_on "fftw"

  def install
    system "./waf", "configure", "--prefix=#{prefix}", "--enable-fftw3"
    system "./waf", "build"
    system "./waf", "install"
  end
end
