VAGRANTFILE_API_VERSION = "2"

Vagrant.configure(VAGRANTFILE_API_VERSION) do |config|
  config.vm.box = "chef/ubuntu-14.04"
  config.vm.provision :shell, path: "bootstrap.sh"
  config.vm.hostname = "tnverify"
  config.vm.network :forwarded_port, guest: 22, host: 11111
  config.vm.synced_folder "/casa3/nobackup/med-m-t/hiseq2/contralateral", "/vagrant/bams"
end
