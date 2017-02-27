# P検定実装

## Usage

```ruby

class Main
  include PValue
  def self.run
    puts PValue.new(Group.new(10,100), Group.new(10, 200)).calculate
  end
end

```

the return value is P_Value.

