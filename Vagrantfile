VAGRANTFILE_API_VERSION = "2"

Vagrant.configure(VAGRANTFILE_API_VERSION) do |config|
  config.vm.box = "chef/ubuntu-14.04"
  config.vm.provision :shell, path: "bootstrap.sh"
  config.name = "tnverify"
  config.cpus = 2
  config.memory = 4096
end
