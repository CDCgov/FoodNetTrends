nextflow.enable.dsl = 2

process TestInputs {
    input:
      val a
      val b
      val c
      val d
      val e
      val f
      val g
      val h
      val i
      val j

    script:
    """
    echo "a: $a"
    echo "b: $b"
    echo "c: $c"
    echo "d: $d"
    echo "e: $e"
    echo "f: $f"
    echo "g: $g"
    echo "h: $h"
    echo "i: $i"
    echo "j: $j"
    """
}

workflow {
  // Pass positional channels without naming them
  TestInputs(
    Channel.value("A"),
    Channel.value("B"),
    Channel.value("C"),
    Channel.value("D"),
    Channel.value("E"),
    Channel.value("F"),
    Channel.value("G"),
    Channel.value("H"),
    Channel.value("I"),
    Channel.value("J")
  )
}

